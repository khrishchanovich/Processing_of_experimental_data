import math
from scipy.stats import chi2, norm
import numpy as np
import matplotlib.pyplot as plt

# В результате эксперимента получены данные, записанные в виде статистического ряда.
data = [0.026, 0.034, 0.028, 0.036, 0.030, 0.038, 0.041, 0.038, 0.030, 0.028,
        0.028, 0.030, 0.034, 0.038, 0.040, 0.036, 0.034, 0.023, 0.032, 0.026,
        0.034, 0.032, 0.024, 0.036, 0.032, 0.026, 0.030, 0.028, 0.038, 0.034,
        0.038, 0.041, 0.028, 0.026, 0.030, 0.034, 0.032, 0.040, 0.036, 0.032,
        0.030, 0.036, 0.034, 0.032, 0.023, 0.032, 0.028, 0.032, 0.026, 0.038,
        0.026, 0.032, 0.028, 0.040, 0.038, 0.030, 0.032, 0.024, 0.036, 0.030,
        0.024, 0.032, 0.030, 0.036, 0.028, 0.041, 0.032, 0.038, 0.034, 0.026,
        0.041, 0.034, 0.023, 0.038, 0.026, 0.030, 0.038, 0.036, 0.040, 0.028,
        0.030, 0.026, 0.034, 0.028, 0.024, 0.036, 0.032, 0.030, 0.038, 0.034,
        0.028, 0.034, 0.040, 0.036, 0.030, 0.038, 0.023, 0.034, 0.032, 0.026
        ]

# 1) Записать значения результатов эксперимента в виде вариационного ряда
data_sort = sorted(data)
print('data_sort', data_sort)

# 2) Найти размах варьирования и разбить его на 9 интервалов
w = round(max(data_sort) - min(data_sort), 3)
print('w', w)

N = len(data_sort)
print('n', N)

k = 9
h = (max(data_sort) - min(data_sort)) / k

print('h', h)

intervals = []
start = min(data_sort)
end = float('{:.3f}'.format(start + h))

while len(intervals) < k:
    interval = [value for value in data_sort if start <= value <= end]
    intervals.append(interval)

    start = end
    end += h
    end = float('{:.3f}'.format(end))

num = 0
for i in intervals:
    num += 1
    if min(i) == max(i):
        print(f'{num}. interval {max(i)}')
    else:
        print(f'{num}. interval {min(i)} - {max(i)}')

print('parts of intervals', intervals)

middle_of_intevals = []
for i in intervals:
    max_val = max(i)
    min_val = min(i)
    middle_of_intevals.append(round((max_val+min_val)/2, 3))

print('middles of intervals', middle_of_intevals)

n = []
count = 0

for i in intervals:
    count = len(i)
    n.append(count)

print('n', n)

Pn = []
for i in n:
    Pn.append(i/h)

print('Pn', Pn)

W = []

for i in n:
    W.append(i/N)

print('W', W)

Pw = []

for i in W:
    Pw.append(i/h)

print('Pw', Pw)

Fx = []
prev_idx = 0

for i in range(len(n)):
    Fx.append((prev_idx + n[i]) / 100)
    prev_idx += n[i]

print('Fx', Fx)

# 3) Построить полигон частот, гистрограмму относительных частот и график эмпирической функции распределения
# plt.figure(figsize=(10, 6))
# plt.plot(middle_of_intevals, n, marker='o', linestyle='-', color='b')
# plt.title('Полигон частот')
# plt.xlabel('Значение интервала')
# plt.ylabel('Частота')
# plt.grid(True)
# plt.show()
#
# plt.figure(figsize=(10, 6))
# plt.bar(middle_of_intevals, Pn, width=h, color='g', edgecolor='black')
# plt.title('Гистограмма относительных частот')
# plt.xlabel('Значение интервала')
# plt.ylabel('Относительная частота')
# plt.grid(True)
# plt.show()
#
# plt.figure(figsize=(10, 6))
# plt.step(middle_of_intevals, Fx, where='post')
# plt.title('График эмпирической функции распределения')
# plt.xlabel('Значение')
# plt.ylabel('Эмпирическая функция распределения')
# plt.grid(True)
# plt.show()

# 4) Найти числовые характеристики выборки xв и Dв
sum = 0
for xi, ni in zip(middle_of_intevals, n):
    sum += xi*ni

x_mid = (1/N) * sum
print('Выборочное среднее', x_mid)

sum = 0
for xi, ni in zip(middle_of_intevals, n):
    sum += (xi**2)*ni

D = (1/N) * sum - x_mid**2
print('Выборочная дисперсия', round(D, 5))

sigma = math.sqrt(D)
print('Стандартное отклонение', sigma)

# 5) Приняв в качестве нулевой гипотeзу H0: генеральная совокупность, из которой
#    извлечена выборка, имеет нормальное распределение, проверть ее, пользуясь
#    критерием Пирсона при уровне значимости alpha = 0.25
z = []
for i in middle_of_intevals:
    z.append(round((i - x_mid)/sigma, 4))
print('z', z)

Fz = []
for zi in z:
    Fz.append(round((1/math.sqrt(2*math.pi)) * math.exp(-(zi**2/2)), 4))
print('Fz', Fz)

temp_n = []
for fzi in Fz:
    temp_n.append(round(((h*N)/sigma) * fzi, 4))
print('temp_n', temp_n)

chi_obs = 0
for ni, temp_ni in zip(n, temp_n):
    chi_obs += ((ni-temp_ni)**2)/temp_ni

alpha = 0.25
df = k - 1
chi_critical = chi2.ppf(1-alpha, df)

print('Наблюдаемое значение', chi_obs)
print('Критическое значение', chi_critical)

if chi_obs < chi_critical:
    print(f'Нет оснований отвергать нулевую гипотезу на уровне значимости {alpha}')
else:
    print(f'Нулевую гипотезу отвергают на уровне значимости {alpha}')

# 6) Найти доверительные интервалы для математического ожидания
#    и среднего квадратичного отклонения при надежности gamma = 0.9
gamma = 0.9
phi_t = gamma / 2

t = 1.65
print('Коэффициент доверия t: ', t)

left_int = round(x_mid - (sigma/math.sqrt(N)) * t, 4)
right_int = round(x_mid + (sigma/math.sqrt(N)) * t, 4)
print(f'Доверительный интервал для математического ожидания: ({left_int}; {right_int})')

q = 0.143
print('Коэффициент q: ', q)

left_int = round(sigma * (1 - q), 4)
right_int = round(sigma * (1 + q), 4)
print(f'Доверительный интервал для среднего квадратичного отклонения ({left_int}; {right_int})')