#include <Wire.h>
#include <math.h>

#define LSM6DS3_ADDR 0x6B

// Регистры LSM6DS3
#define WHO_AM_I   0x0F
#define CTRL1_XL   0x10  // Акселерометр
#define CTRL2_G    0x11  // Гироскоп
#define CTRL3_C    0x12  // Управление
#define STATUS_REG 0x1E  // Статус готовности данных
#define OUTX_L_G   0x22  // Выходные данные гироскопа и акселерометра

// Коэффициент чувствительности для шкалы ±500 dps (17.50 mdps/LSB) в рад/с:
const float GYRO_SCALE_RAD = 0.000305432619f;
// Коэффициент чувствительности для акселерометра ±4g (0.122 mg/LSB) в g:
const float ACCEL_SCALE_G  = 0.000122f;

// Смещения нуля гироскопа
float gyroBiasX = 0.0f;
float gyroBiasY = 0.0f;
float gyroBiasZ = 0.0f;

// Масштаб вектора гравитации (калибруется в покое)
float gravityMag = 1.0f;

// Смещения нуля линейного ускорения (в см/с^2)
float accelBiasX = 0.0f;
float accelBiasY = 0.0f;
float accelBiasZ = 0.0f;

// Позиция (в сантиметрах) и линейная скорость (в см/с)
float posX = 0.0f, posY = 0.0f, posZ = 0.0f;
float velX = 0.0f, velY = 0.0f, velZ = 0.0f;

// Счетчик стабильного покоя для алгоритма ZUPT
int restCounter = 0;

// Кватернион ориентации: Q = (q0, q1, q2, q3)
float q0 = 1.0f, q1 = 0.0f, q2 = 0.0f, q3 = 0.0f;

// Коэффициент пропорциональной коррекции горизонта Mahony
const float Kp = 0.5f;

// Фильтр низких частот (PT1 LPF) для гироскопа (срез 25 Гц)
float lpf_gx = 0.0f, lpf_gy = 0.0f, lpf_gz = 0.0f;
const float LPF_CUTOFF_HZ = 25.0f;

unsigned long prevTime = 0;
unsigned long lastTelemetryTime = 0;

// Последние сырые показания для телеметрии
float last_rax = 0.0f, last_ray = 0.0f, last_raz = 0.0f;
float last_rgx = 0.0f, last_rgy = 0.0f, last_rgz = 0.0f;

void writeReg(uint8_t reg, uint8_t val) {
  Wire.beginTransmission(LSM6DS3_ADDR);
  Wire.write(reg);
  Wire.write(val);
  Wire.endTransmission();
}

// Прямое чтение 12 байт датчика без искажений
bool readSensorsRawData(int16_t &gX, int16_t &gY, int16_t &gZ, int16_t &aX, int16_t &aY, int16_t &aZ) {
  Wire.beginTransmission(LSM6DS3_ADDR);
  Wire.write(OUTX_L_G);
  if (Wire.endTransmission(false) != 0) return false;

  if (Wire.requestFrom((uint8_t)LSM6DS3_ADDR, (uint8_t)12) == 12) {
    uint8_t b[12];
    for (int i = 0; i < 12; i++) {
      b[i] = Wire.read();
    }
    gX = (int16_t)(b[0] | (b[1] << 8));
    gY = (int16_t)(b[2] | (b[3] << 8));
    gZ = (int16_t)(b[4] | (b[5] << 8));
    aX = (int16_t)(b[6] | (b[7] << 8));
    aY = (int16_t)(b[8] | (b[9] << 8));
    aZ = (int16_t)(b[10] | (b[11] << 8));
    return true;
  }
  return false;
}

// Прецизионная калибровка нуля гироскопа и вектора гравитации
void calibrate() {
  Serial.println(F("[CALIB] Калибровка гироскопа и акселерометра... Держите плату неподвижно!"));

  bool calibrated = false;
  while (!calibrated) {
    float sumX = 0.0f, sumY = 0.0f, sumZ = 0.0f;
    float sumAx = 0.0f, sumAy = 0.0f, sumAz = 0.0f;
    float minZ = 1e9f, maxZ = -1e9f;
    const int SAMPLES = 400;

    for (int i = 0; i < SAMPLES; i++) {
      int16_t gX, gY, gZ, aX, aY, aZ;
      if (readSensorsRawData(gX, gY, gZ, aX, aY, aZ)) {
        float gx = gX * GYRO_SCALE_RAD;
        float gy = gY * GYRO_SCALE_RAD;
        float gz = gZ * GYRO_SCALE_RAD;

        sumX += gx;
        sumY += gy;
        sumZ += gz;

        sumAx += aX * ACCEL_SCALE_G;
        sumAy += aY * ACCEL_SCALE_G;
        sumAz += aZ * ACCEL_SCALE_G;

        if (gz < minZ) minZ = gz;
        if (gz > maxZ) maxZ = gz;
      }
      delay(2);
    }

    // Допуск по шуму Z: 0.08 рад/с (~4.6 dps)
    if ((maxZ - minZ) < 0.08f) {
      gyroBiasX = sumX / SAMPLES;
      gyroBiasY = sumY / SAMPLES;
      gyroBiasZ = sumZ / SAMPLES;

      float meanAx = sumAx / SAMPLES;
      float meanAy = sumAy / SAMPLES;
      float meanAz = sumAz / SAMPLES;
      gravityMag = sqrtf(meanAx * meanAx + meanAy * meanAy + meanAz * meanAz);
      if (gravityMag < 0.8f || gravityMag > 1.2f) gravityMag = 1.0f;

      calibrated = true;
      Serial.printf("[CALIB] Успех! Bias (dps): X=%.2f, Y=%.2f, Z=%.2f | G=%.3fg\n",
                    gyroBiasX * (180.0f / M_PI),
                    gyroBiasY * (180.0f / M_PI),
                    gyroBiasZ * (180.0f / M_PI),
                    gravityMag);
    } else {
      Serial.println(F("[CALIB] Движение обнаружено! Повтор через 100 мс..."));
      delay(100);
    }
  }

  q0 = 1.0f; q1 = 0.0f; q2 = 0.0f; q3 = 0.0f;
  lpf_gx = 0.0f; lpf_gy = 0.0f; lpf_gz = 0.0f;

  posX = 0.0f; posY = 0.0f; posZ = 0.0f;
  velX = 0.0f; velY = 0.0f; velZ = 0.0f;
  accelBiasX = 0.0f; accelBiasY = 0.0f; accelBiasZ = 0.0f;
  restCounter = 20;
}

// 6-осевой фильтр Mahony AHRS с динамической изоляцией курса
void mahonyAHRSUpdate(float gx, float gy, float gz, float ax, float ay, float az, float dt) {
  float aMag = sqrtf(ax * ax + ay * ay + az * az);

  if (aMag > 0.05f) {
    float recipMag = 1.0f / aMag;
    float norm_ax = ax * recipMag;
    float norm_ay = ay * recipMag;
    float norm_az = az * recipMag;

    float vx = 2.0f * (q1 * q3 - q0 * q2);
    float vy = 2.0f * (q0 * q1 + q2 * q3);
    float vz = q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3;

    float ex = norm_ay * vz - norm_az * vy;
    float ey = norm_az * vx - norm_ax * vz;
    float ez = norm_ax * vy - norm_ay * vx;

    float aDev = fabsf(aMag - 1.0f);
    float accelWeight = 1.0f;
    if (aDev > 0.10f) {
      if (aDev >= 0.25f) {
        accelWeight = 0.0f;
      } else {
        accelWeight = (0.25f - aDev) / 0.15f;
      }
    }

    gx += (Kp * accelWeight) * ex;
    gy += (Kp * accelWeight) * ey;
    gz += (Kp * accelWeight) * ez;
  }

  float qDot1 = 0.5f * (-q1 * gx - q2 * gy - q3 * gz);
  float qDot2 = 0.5f * ( q0 * gx + q2 * gz - q3 * gy);
  float qDot3 = 0.5f * ( q0 * gy - q1 * gz + q3 * gx);
  float qDot4 = 0.5f * ( q0 * gz + q1 * gy - q2 * gx);

  q0 += qDot1 * dt;
  q1 += qDot2 * dt;
  q2 += qDot3 * dt;
  q3 += qDot4 * dt;

  float recipNorm = 1.0f / sqrtf(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
  q0 *= recipNorm;
  q1 *= recipNorm;
  q2 *= recipNorm;
  q3 *= recipNorm;
}

// Преобразование кватерниона в углы Эйлера (Pitch, Roll, Yaw)
void getEuler(float &roll, float &pitch, float &yaw) {
  float sin_a = 2.0f * (q0 * q1 + q2 * q3);
  float cos_a = 1.0f - 2.0f * (q1 * q1 + q2 * q2);
  pitch = atan2f(sin_a, cos_a) * (180.0f / M_PI);

  float sin_b = 2.0f * (q0 * q2 - q3 * q1);
  if (fabsf(sin_b) >= 1.0f)
    roll = copysignf(90.0f, sin_b);
  else
    roll = asinf(sin_b) * (180.0f / M_PI);

  float sin_c = 2.0f * (q0 * q3 + q1 * q2);
  float cos_c = 1.0f - 2.0f * (q2 * q2 + q3 * q3);
  yaw = atan2f(sin_c, cos_c) * (180.0f / M_PI);
}

void setup() {
  Serial.begin(115200);
  while (!Serial && millis() < 2500);

  Wire.begin(8, 9);
  Wire.setClock(400000);
  delay(100);

  // BDU = 1, IF_INC = 1
  writeReg(CTRL3_C, 0x44);
  // Акселерометр: 416 Гц, ±4g
  writeReg(CTRL1_XL, 0x68);
  // Гироскоп: 416 Гц, ±500 dps
  writeReg(CTRL2_G, 0x64);

  delay(50);
  calibrate();

  prevTime = micros();
  lastTelemetryTime = millis();
}

void loop() {
  if (Serial.available()) {
    char cmd = Serial.read();
    if (cmd == 't' || cmd == 'T') {
      q0 = 1.0f; q1 = 0.0f; q2 = 0.0f; q3 = 0.0f;
    } else if (cmd == 'c' || cmd == 'C') {
      calibrate();
    } else if (cmd == 'r' || cmd == 'R' || cmd == '0') {
      // Сброс координат позиции и скорости в ноль
      posX = 0.0f; posY = 0.0f; posZ = 0.0f;
      velX = 0.0f; velY = 0.0f; velZ = 0.0f;
    }
  }

  unsigned long now = micros();
  float dt = (now - prevTime) / 1000000.0f;
  prevTime = now;
  if (dt <= 0.0f || dt > 0.05f) dt = 0.0025f;

  int16_t gX, gY, gZ, aX, aY, aZ;
  if (readSensorsRawData(gX, gY, gZ, aX, aY, aZ)) {
    // Вычитание смещения гироскопа (в рад/с)
    float raw_gx = (gX * GYRO_SCALE_RAD) - gyroBiasX;
    float raw_gy = (gY * GYRO_SCALE_RAD) - gyroBiasY;
    float raw_gz = (gZ * GYRO_SCALE_RAD) - gyroBiasZ;

    // Акселерометр в g
    float raw_ax = aX * ACCEL_SCALE_G;
    float raw_ay = aY * ACCEL_SCALE_G;
    float raw_az = aZ * ACCEL_SCALE_G;

    // LPF фильтр 25 Гц
    float alpha = (2.0f * M_PI * LPF_CUTOFF_HZ * dt) / (1.0f + 2.0f * M_PI * LPF_CUTOFF_HZ * dt);
    lpf_gx += alpha * (raw_gx - lpf_gx);
    lpf_gy += alpha * (raw_gy - lpf_gy);
    lpf_gz += alpha * (raw_gz - lpf_gz);

    // Согласование осей фильтра:
    // a1 = Вперед (-raw_ay), a2 = Вправо (raw_ax), a3 = Вверх (raw_az)
    float a1 = -raw_ay;
    float a2 =  raw_ax;
    float a3 =  raw_az;

    float g1 = -lpf_gy;
    float g2 =  lpf_gx;
    float g3 =  lpf_gz;

    // -------------------------------------------------------------
    // ИНЕРЦИАЛЬНАЯ НАВИГАЦИЯ (SINS): ДВОЙНОЕ ИНТЕГРИРОВАНИЕ В МИРЕ
    // -------------------------------------------------------------
    // 1. Поворот вектора ускорения из корпуса в мировую систему по кватерниону Q:
    float q0q0 = q0 * q0;
    float q1q1 = q1 * q1;
    float q2q2 = q2 * q2;
    float q3q3 = q3 * q3;

    float r11 = q0q0 + q1q1 - q2q2 - q3q3;
    float r12 = 2.0f * (q1 * q2 - q0 * q3);
    float r13 = 2.0f * (q1 * q3 + q0 * q2);

    float r21 = 2.0f * (q1 * q2 + q0 * q3);
    float r22 = q0q0 - q1q1 + q2q2 - q3q3;
    float r23 = 2.0f * (q2 * q3 - q0 * q1);

    float r31 = 2.0f * (q1 * q3 - q0 * q2);
    float r32 = 2.0f * (q2 * q3 + q0 * q1);
    float r33 = q0q0 - q1q1 - q2q2 + q3q3;

    float a_world_x = r11 * a1 + r12 * a2 + r13 * a3;
    float a_world_y = r21 * a1 + r22 * a2 + r23 * a3;
    float a_world_z = r31 * a1 + r32 * a2 + r33 * a3;

    // 2. Вычитание мировой гравитации (по оси Z вверх)
    float a_lin_x = a_world_x;
    float a_lin_y = a_world_y;
    float a_lin_z = a_world_z - gravityMag;

    // 3. Перевод в см/с^2 (1g = 980.665 см/с^2)
    const float G_TO_CMS2 = 980.665f;
    float linAccX = a_lin_x * G_TO_CMS2 - accelBiasX;
    float linAccY = a_lin_y * G_TO_CMS2 - accelBiasY;
    float linAccZ = a_lin_z * G_TO_CMS2 - accelBiasZ;

    // Мертвая зона по шуму АЦП (~6-8 см/с^2)
    if (fabsf(linAccX) < 7.0f) linAccX = 0.0f;
    if (fabsf(linAccY) < 7.0f) linAccY = 0.0f;
    if (fabsf(linAccZ) < 8.5f) linAccZ = 0.0f;

    // 4. Детектор покоя (Zero Velocity Detection - ZUPT)
    float gyroSpeed = sqrtf(g1 * g1 + g2 * g2 + g3 * g3);
    float linAccMag = sqrtf(linAccX * linAccX + linAccY * linAccY + linAccZ * linAccZ);

    bool isRest = (gyroSpeed < 0.025f) && (linAccMag < 14.0f);
    if (isRest) {
      if (restCounter < 25) restCounter++;
    } else {
      restCounter = 0;
    }

    if (restCounter >= 10) {
      // ПОКОЙ: скорость строго обнуляется, устраняя экспоненциальный дрейф позиции
      velX = 0.0f;
      velY = 0.0f;
      velZ = 0.0f;

      // Мягкая подстройка смещения нуля гироскопа и акселерометра в покое
      gyroBiasX += raw_gx * 0.0003f;
      gyroBiasY += raw_gy * 0.0003f;
      gyroBiasZ += raw_gz * 0.0003f;

      accelBiasX += a_lin_x * G_TO_CMS2 * 0.002f;
      accelBiasY += a_lin_y * G_TO_CMS2 * 0.002f;
      accelBiasZ += a_lin_z * G_TO_CMS2 * 0.002f;

      g1 = 0.0f; g2 = 0.0f; g3 = 0.0f;
    } else {
      // ДВИЖЕНИЕ: интегрируем ускорение в скорость и перемещение
      velX += linAccX * dt;
      velY += linAccY * dt;
      velZ += linAccZ * dt;

      // Мягкое экспоненциальное затухание скорости (damping), стабилизирующее свободный полет
      velX *= (1.0f - 0.35f * dt);
      velY *= (1.0f - 0.35f * dt);
      velZ *= (1.0f - 0.35f * dt);

      posX += velX * dt;
      posY += velY * dt;
      posZ += velZ * dt;
    }

    // Обновление ориентации Mahony AHRS
    mahonyAHRSUpdate(g1, g2, g3, a1, a2, a3, dt);

    last_rax = raw_ax;
    last_ray = raw_ay;
    last_raz = raw_az;
    last_rgx = lpf_gx;
    last_rgy = lpf_gy;
    last_rgz = lpf_gz;
  }

  // Передача пакета телеметрии каждые 30 мс
  if (millis() - lastTelemetryTime >= 30) {
    lastTelemetryTime = millis();

    float roll, pitch, yaw;
    getEuler(roll, pitch, yaw);

    float gx_dps = last_rgx * (180.0f / M_PI);
    float gy_dps = last_rgy * (180.0f / M_PI);
    float gz_dps = last_rgz * (180.0f / M_PI);

    Serial.printf("Roll:%.2f,Pitch:%.2f,Yaw:%.2f,AX:%.2f,AY:%.2f,AZ:%.2f,GX:%.1f,GY:%.1f,GZ:%.1f,PX:%.1f,PY:%.1f,PZ:%.1f,VX:%.1f,VY:%.1f,VZ:%.1f,Q:%.4f,%.4f,%.4f,%.4f\n",
                  roll, pitch, yaw, last_rax, last_ray, last_raz, gx_dps, gy_dps, gz_dps,
                  posX, posY, posZ, velX, velY, velZ, q0, q1, q2, q3);
  }

  delayMicroseconds(2400);
}
