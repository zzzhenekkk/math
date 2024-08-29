#include "s21_math.h"

// возвращает ближайшее целое число, не меньшее заданного значения по модулю
long double s21_ceil(double x) {
  long long int res = s21_floor(x) + (x != (long long int)x);
  return (long double)((x > 1e17 || x < -1e17 || s21_isnan(x) || x == -0.)
                           ? x
                           : res);
}

// вычисляет абсолютное значение целого числа
int s21_abs(int x) {
  if (x < 0) x = -x;
  return x;
}

// вычисляет арксинус
long double s21_asin(double x) {
  long double res = 0.0;
  if (s21_fabs(x) == 1)
    res = x * PI / 2;
  else if (s21_fabs(x) < 1) {
    double n = 0.0000001, asn = 1.0;
    int i = 0;
    while (s21_fabs(asn) > n) {
      asn = s21_factorial(2 * i) * s21_pow(x, 2 * i + 1) /
            (s21_pow(4, i) * s21_pow(s21_factorial(i), 2) * (2 * i + 1));
      res += asn;
      i++;
    }
  } else
    res = S21_NAN;
  return res;
}

// вычисляет арккосинус
long double s21_acos(double x) {
  long double res = S21_NAN;
  if (s21_fabs(x) < 1)
    res = PI / 2 - s21_asin(x);
  else if (x == 1)
    res = 0;
  else if (x == -1)
    res = PI;
  return res;
}

// вычисляет арктангенс
long double s21_atan(double x) {
  long double res = S21_NAN;
  if (x == S21_INF) return PI / 2;
  if (x == -S21_INF) return -PI / 2;
  if (x != x) return res;
  double min = -PI / 2, max = PI / 2;
  while (min < max) {
    res = (max - min) / 2 + min;
    long double tn = s21_tan(res);
    if (tn == x)
      min = max;
    else if (tn < x)
      min = res + 0.00000001;
    else
      max = res - 0.00000001;
  }
  return res;
}

// вычисляет косинус
long double s21_cos(double x) { return s21_sin(PI / 2 + x); }

// возвращает значение e, возведенное в заданную степень
long double s21_exp(double x) {
  long double prev = 0., cur = 1., y = (long double)x;
  long double k = 1., den = 1.;
  int flag = 0;
  if (x < 0) {
    x = -x;
    y = -y;
    flag = 1;
  }
  while (s21_fabs(cur - prev) > EPS) {
    prev = cur;
    cur += y / den;
    y *= x;
    den *= (k + 1);
    k++;
  }
  return (flag) ? 1 / cur : cur;
}

// вычисляет абсолютное значение числа с плавающей точкой
long double s21_fabs(double x) {
  long double res = x;
  if (x == 0)
    res = 0;
  else if (x < 0)
    res = -x;
  return res;
}

// возвращает ближайшее целое число, не превышающее заданное значение
long double s21_floor(double x) {
  long long int res = (long long int)x;
  if (x < 0 && x != (double)res) res--;
  return (long double)((x > 1e17 || x < -1e17 || s21_isnan(x) || x == -0.)
                           ? x
                           : res);
}

// остаток операции деления с плавающей точкой
long double s21_fmod(double x, double y) {
  long double x_mod = s21_fabs(x);
  long double y_mod = s21_fabs(y);
  long double dif = x_mod;
  if (y_mod == 0 || s21_isinf(x) || s21_isnan(y_mod) || s21_isnan(x_mod)) {
    dif = S21_NAN;
  } else if (s21_isinf(y_mod))
    dif = x;
  else {
    dif = x_mod / y_mod - s21_floor(x_mod / y_mod);
    dif *= y_mod;
    if (x < 0) dif = -dif;
  }
  return dif;
}

// возводит число в указанную степень
long double s21_pow(double base, double exp) {
  long double res = 1;
  int sign = 1;
  if (exp == 0)
    res = 1;
  else if (s21_isnan(exp))
    res = S21_NAN;
  else if (base == 0) {
    if (exp == -S21_INF || exp < 0)
      res = S21_INF;
    else
      res = 0;
  } else if (s21_isnan(base))
    res = S21_NAN;
  else if (s21_isinf(base)) {
    if (exp == -S21_INF || exp < 0)
      res = 0;
    else
      res = S21_INF;
  } else if (s21_isinf(exp)) {
    if (exp == S21_INF)
      res = S21_INF;
    else
      res = 0;
  } else if (base < 0 && s21_floor(exp) != exp)
    res = S21_NAN;
  else {
    if (base < 0) sign *= -1;
    res = s21_exp(exp * s21_log(base * sign));
    if (s21_fmod(exp, 2)) res *= sign;
  }
  return res;
}

// вычисляет квадратный корень из числа
long double s21_sqrt(double x) {
  long double res;
  if (s21_isnan(x) || x == -S21_INF || x < 0)
    res = S21_NAN;
  else if (s21_isinf(x))
    res = S21_INF;
  else
    res = s21_pow(x, 0.5);
  return res;
}

// вычисляет логарифм
long double s21_log(double x) {
  int e_count = 0;
  double start = x;
  long double res = 0.0, q = -1;
  if (x < 0 || s21_isnan(x))
    res = S21_NAN;
  else if (x == 0)
    res = S21_INF;
  else if (s21_isinf(x))
    res = S21_INF;
  else {
    if (0 < x && x < 1) {
      while (x <= s21_e) {
        x *= s21_e;
        e_count++;
      }
    } else {
      while (x >= s21_e) {
        x /= s21_e;
        e_count++;
      }
    }
    while (s21_exp(res) != x) {
      if (res - q < 0.0000001) break;
      q = res;
      res = q + 2 * (x - s21_exp(q)) / (x + s21_exp(q));
    }
    if (0 < start && start < 1)
      res -= e_count;
    else
      res += e_count;
  }
  return res;
}

// вычисляет синус
long double s21_sin(double x) {
  long double s = 0;
  double q = 1, n = 0.000001;
  int i = 0;
  if (s21_isinf(x) || s21_isnan(x))
    s = S21_NAN;
  else {
    x = s21_fmod(x, 2 * PI);
    while (s21_fabs(q) >= n) {
      q = s21_pow(-1, i) * s21_pow(x, 2 * i + 1) / s21_factorial(2 * i + 1);
      s += q;
      i++;
    }
  }
  return s;
}

// вычисляет тангенс
long double s21_tan(double x) {
  long double res = 0;
  res = s21_sin(x) / s21_cos(x);
  return res;
}

// вспомогательные функции

long double s21_factorial(double x) {
  double s = 1.0;
  if (x >= 0)
    for (int i = 1; i <= x; i++) s *= i;
  else
    s = 0;
  return s;
}

int s21_isnan(long double x) { return (x != x) ? 1 : 0; }
// проверяет число на nan

int s21_isinf(long double x) {
  return ((x == S21_INF || x == -S21_INF) && !s21_isnan(x)) ? 1 : 0;
}
// проверяет число на inf

long double s21_tolerance(long double x) {
  while (x >= 1e16 && x < S21_INF) x /= 10;
  return x;
}
