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
// 0.01750 * (PI / 180) = 0.000305432619f
const float GYRO_SCALE_RAD = 0.000305432619f;
// Коэффициент чувствительности для акселерометра ±4g (0.122 mg/LSB) в g:
const float ACCEL_SCALE_G  = 0.000122f;

// Смещения нуля гироскопа
float gyroBiasX = 0.0f;
float gyroBiasY = 0.0f;
float gyroBiasZ = 0.0f;

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

// Прецизионная калибровка нуля гироскопа при старте с контролем покоя
void calibrate() {
  Serial.println(F("[CALIB] Калибровка гироскопа... Держите плату неподвижно!"));

  bool calibrated = false;
  while (!calibrated) {
    float sumX = 0.0f, sumY = 0.0f, sumZ = 0.0f;
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

        if (gz < minZ) minZ = gz;
        if (gz > maxZ) maxZ = gz;
      }
      delay(2);
    }

    // Допуск по шуму Z: 0.08 рад/с (~4.6 dps) - мягкий порог для устойчивой калибровки
    if ((maxZ - minZ) < 0.08f) {
      gyroBiasX = sumX / SAMPLES;
      gyroBiasY = sumY / SAMPLES;
      gyroBiasZ = sumZ / SAMPLES;
      calibrated = true;
      Serial.printf("[CALIB] Успех! Bias (dps): X=%.2f, Y=%.2f, Z=%.2f\n",
                    gyroBiasX * (180.0f / M_PI),
                    gyroBiasY * (180.0f / M_PI),
                    gyroBiasZ * (180.0f / M_PI));
    } else {
      Serial.println(F("[CALIB] Движение обнаружено! Повтор через 100 мс..."));
      delay(100);
    }
  }

  q0 = 1.0f; q1 = 0.0f; q2 = 0.0f; q3 = 0.0f;
  lpf_gx = 0.0f; lpf_gy = 0.0f; lpf_gz = 0.0f;
}

// 6-осевой фильтр Mahony AHRS с динамической изоляцией курса при наклонах и буграх
void mahonyAHRSUpdate(float gx, float gy, float gz, float ax, float ay, float az, float dt) {
  // Вычисляем модуль вектора ускорения
  float aMag = sqrtf(ax * ax + ay * ay + az * az);

  // Коррекция горизонта по акселерометру
  if (aMag > 0.05f) {
    // Нормализация вектора акселерометра
    float recipMag = 1.0f / aMag;
    float norm_ax = ax * recipMag;
    float norm_ay = ay * recipMag;
    float norm_az = az * recipMag;

    // Вектор гравитации в системе координат датчика по текущему кватерниону
    float vx = 2.0f * (q1 * q3 - q0 * q2);
    float vy = 2.0f * (q0 * q1 + q2 * q3);
    float vz = q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3;

    // Ошибка гравитации: векторное произведение измеренного ускорения и расчетной вертикали:
    // e = a_norm x v
    // Так как e ортогонален v (вертикали), поправка по оси курса (Yaw) тождественно равна 0!
    float ex = norm_ay * vz - norm_az * vy;
    float ey = norm_az * vx - norm_ax * vz;
    float ez = norm_ax * vy - norm_ay * vx;

    // ДИНАМИЧЕСКОЕ ПОДАВЛЕНИЕ УСКОРЕНИЙ (Acceleration Rejection):
    // При езде по буграм, толчках или разгоне робота aMag отклоняется от 1.0g.
    // Если модуль ускорения отличается от 1g более чем на 0.15g, снижаем вес акселерометра,
    // чтобы толчки шасси не наклоняли виртуальный горизонт.
    float aDev = fabsf(aMag - 1.0f);
    float accelWeight = 1.0f;
    if (aDev > 0.10f) {
      if (aDev >= 0.25f) {
        accelWeight = 0.0f; // Полный отказ от акселерометра во время резкого толчка/кочки
      } else {
        accelWeight = (0.25f - aDev) / 0.15f; // Плавное снижение доверия
      }
    }

    // Применение пропорциональной коррекции к гироскопу
    gx += (Kp * accelWeight) * ex;
    gy += (Kp * accelWeight) * ey;
    gz += (Kp * accelWeight) * ez;
  }

  // Численное интегрирование кинематического уравнения кватерниона:
  // dQ/dt = 0.5 * Q (x) omega
  float qDot1 = 0.5f * (-q1 * gx - q2 * gy - q3 * gz);
  float qDot2 = 0.5f * ( q0 * gx + q2 * gz - q3 * gy);
  float qDot3 = 0.5f * ( q0 * gy - q1 * gz + q3 * gx);
  float qDot4 = 0.5f * ( q0 * gz + q1 * gy - q2 * gx);

  q0 += qDot1 * dt;
  q1 += qDot2 * dt;
  q2 += qDot3 * dt;
  q3 += qDot4 * dt;

  // Нормализация кватерниона
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

  // Настройка LSM6DS3:
  // CTRL3_C: BDU = 1 (бит 6, защита от склеивания байтов), IF_INC = 1 (бит 2, автоинкремент)
  writeReg(CTRL3_C, 0x44);

  // CTRL1_XL: 416 Гц ODR, диапазон ±4g, фильтр 400 Гц
  writeReg(CTRL1_XL, 0x68);

  // CTRL2_G: 416 Гц ODR, диапазон ±500 dps (FS=01) для максимальной точности на роботе
  writeReg(CTRL2_G, 0x64);

  delay(50);
  calibrate();

  prevTime = micros();
  lastTelemetryTime = millis();
}

void loop() {
  // Обработка команд отладки и визуализатора по Serial
  if (Serial.available()) {
    char cmd = Serial.read();
    if (cmd == 't' || cmd == 'T' || cmd == '0') {
      q0 = 1.0f; q1 = 0.0f; q2 = 0.0f; q3 = 0.0f;
    } else if (cmd == 'c' || cmd == 'C') {
      calibrate();
    }
  }

  unsigned long now = micros();
  float dt = (now - prevTime) / 1000000.0f;
  prevTime = now;
  if (dt <= 0.0f || dt > 0.05f) dt = 0.0025f;

  int16_t gX, gY, gZ, aX, aY, aZ;
  if (readSensorsRawData(gX, gY, gZ, aX, aY, aZ)) {
    // Перевод в физические единицы:
    // Гироскоп с вычетом откалиброванного смещения (рад/с):
    float raw_gx = (gX * GYRO_SCALE_RAD) - gyroBiasX;
    float raw_gy = (gY * GYRO_SCALE_RAD) - gyroBiasY;
    float raw_gz = (gZ * GYRO_SCALE_RAD) - gyroBiasZ;

    // Акселерометр в g:
    float raw_ax = aX * ACCEL_SCALE_G;
    float raw_ay = aY * ACCEL_SCALE_G;
    float raw_az = aZ * ACCEL_SCALE_G;

    // Фильтр низких частот PT1 (срез 25 Гц) для подавления вибраций шасси
    float alpha = (2.0f * M_PI * LPF_CUTOFF_HZ * dt) / (1.0f + 2.0f * M_PI * LPF_CUTOFF_HZ * dt);
    lpf_gx += alpha * (raw_gx - lpf_gx);
    lpf_gy += alpha * (raw_gy - lpf_gy);
    lpf_gz += alpha * (raw_gz - lpf_gz);

    // -------------------------------------------------------------
    // МЯГКИЙ АНТИ-ДРЕЙФ В ПОКОЕ (Static Rest ZUPT):
    // Применяется только при подтвержденном абсолютном покое на столе:
    // акселерометр стабилен (~1.0g), а отфильтрованная скорость близка к нулю.
    // -------------------------------------------------------------
    float accelMag = sqrtf(raw_ax * raw_ax + raw_ay * raw_ay + raw_az * raw_az);
    float gyroMag  = sqrtf(lpf_gx * lpf_gx + lpf_gy * lpf_gy + lpf_gz * lpf_gz);

    float gz_for_filter = lpf_gz;
    float gx_for_filter = lpf_gx;
    float gy_for_filter = lpf_gy;

    if (fabsf(accelMag - 1.0f) < 0.035f && gyroMag < 0.015f) { // < ~0.85 dps
      // Медленная адаптация смещения нуля под температуру кристалла
      gyroBiasX += raw_gx * 0.0003f;
      gyroBiasY += raw_gy * 0.0003f;
      gyroBiasZ += raw_gz * 0.0003f;

      // Обнуление угловой скорости в покое исключает любой уход курса на столе
      gx_for_filter = 0.0f;
      gy_for_filter = 0.0f;
      gz_for_filter = 0.0f;
    }

    // Передача в фильтр Mahony AHRS с согласованием осей корпуса:
    // Ось 1 (Pitch, наклон носа): -gy
    // Ось 2 (Roll, крен вправо):   gx
    // Ось 3 (Yaw, курс):           gz
    mahonyAHRSUpdate(-gy_for_filter, gx_for_filter, gz_for_filter, -raw_ay, raw_ax, raw_az, dt);

    last_rax = raw_ax;
    last_ray = raw_ay;
    last_raz = raw_az;
    last_rgx = lpf_gx;
    last_rgy = lpf_gy;
    last_rgz = lpf_gz;
  }

  // Передача пакета телеметрии каждые 30 мс (совместимо с visualizer.html)
  if (millis() - lastTelemetryTime >= 30) {
    lastTelemetryTime = millis();

    float roll, pitch, yaw;
    getEuler(roll, pitch, yaw);

    float gx_dps = last_rgx * (180.0f / M_PI);
    float gy_dps = last_rgy * (180.0f / M_PI);
    float gz_dps = last_rgz * (180.0f / M_PI);

    Serial.printf("Roll:%.2f,Pitch:%.2f,Yaw:%.2f,AX:%.2f,AY:%.2f,AZ:%.2f,GX:%.1f,GY:%.1f,GZ:%.1f,Q:%.4f,%.4f,%.4f,%.4f\n",
                  roll, pitch, yaw, last_rax, last_ray, last_raz, gx_dps, gy_dps, gz_dps, q0, q1, q2, q3);
  }

  delayMicroseconds(2400); // Согласовано с частотой ODR 416 Гц
}
