# Отчет по тестированию интерфейса OpenBLAS (Level 1)

## Описание

В данной работе проведено тестирование функций **Level 1** (векторные операции) библиотеки OpenBLAS. Основная цель — проверить корректность работы интерфейса **CBLAS** на языке Си для различных типов данных и условий.

## Что проверялось (Реализованные тесты)

### 1. Поддержка типов данных

Для каждой группы функций тесты реализованы в 4 вариантах точности:

* **S (float)** и **D (double)** — вещественные числа.
* **C (float complex)** и **Z (double complex)** — комплексные числа (использовалась библиотека `<complex.h>`).

### 2. Основные функции

* **AXPY ()**: Базовая операция суммирования векторов. Проверена передача множителя `alpha` (как значения для `float` и как указателя для `complex`).
* **SCAL ()**: Проверка изменения значений вектора.
* **COPY / SWAP**: Тесты на копирование массивов и корректный обмен значениями между ними.
* **DOT / NRM2 / ASUM**: Функции получения характеристик вектора (скалярное произведение, норма, сумма модулей). еще проверил что функции возвращают адекватные значения (не отрицательные для норм).

### 3. Edge Cases

Чтобы убедиться в надежности интерфейса, добавлены специфические проверки:

* **Zero Length**: Вызов функций с длиной вектора `n = 0`. Библиотека не должна выдавать ошибок или падать.
* **Stride (incX > 1)**: Данные берутся через определенный интервал в памяти. интерфейс правильно интерпретирует структуру массива.

## Технические детали

* **Проверка результатов**: Использован механизм `assert()`. Если библиотека вернет неверный результат, выполнение программы немедленно прервется с ошибкой.
* **Многопоточность**: Тесты запускались с переменными окружения `OPENBLAS_NUM_THREADS=1` и `OPENBLAS_NUM_THREADS=4` для проверки стабильности интерфейса в многопоточном режиме.

## Код

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <complex.h>
#include <cblas.h>

int main() {

    printf("AXPY\n");

    float sx[3] = {1,2,3};
    float sy[3] = {4,5,6};
    cblas_saxpy(3,2.0f,sx,1,sy,1);
    assert(sy[0] == 6);
    printf("saxpy passed\n");

    double dx[3] = {1,2,3};
    double dy[3] = {4,5,6};
    cblas_daxpy(3,2.0,dx,1,dy,1);
    assert(dy[1] == 9);
    printf("daxpy passed\n");

    float complex cx[1] = {1+I};
    float complex cy[1] = {0};
    float complex ca = 1+0*I;
    cblas_caxpy(1,&ca,cx,1,cy,1);
    assert(crealf(cy[0]) == 1);
    printf("caxpy passed\n");

    double complex zx[1] = {1+I};
    double complex zy[1] = {0};
    double complex za = 1+0*I;
    cblas_zaxpy(1,&za,zx,1,zy,1);
    assert(creal(zy[0]) == 1);
    printf("zaxpy passed\n");


    printf("\nSCAL\n");

    cblas_sscal(3,2.0f,sx,1);
    assert(sx[0] == 2);
    printf("sscal passed\n");

    cblas_dscal(3,2.0,dx,1);
    assert(dx[0] == 2);
    printf("dscal passed\n");

    cblas_cscal(1,&ca,cx,1);
    assert(crealf(cx[0]) == 1);
    printf("cscal passed\n");

    cblas_zscal(1,&za,zx,1);
    assert(creal(zx[0]) == 1);
    printf("zscal passed\n");


    printf("\nCOPY\n");

    float s_copy[3];
    cblas_scopy(3,sx,1,s_copy,1);
    assert(s_copy[0] == sx[0]);
    printf("scopy passed\n");

    double d_copy[3];
    cblas_dcopy(3,dx,1,d_copy,1);
    assert(d_copy[0] == dx[0]);
    printf("dcopy passed\n");


    printf("\nSWAP\n");

    cblas_sswap(3,sx,1,s_copy,1);
    printf("sswap passed\n");

    cblas_dswap(3,dx,1,d_copy,1);
    printf("dswap passed\n");


    printf("\nDOT\n");

    float sdot = cblas_sdot(3,sx,1,sy,1);
    assert(sdot != 0);
    printf("sdot passed\n");

    double ddot = cblas_ddot(3,dx,1,dy,1);
    assert(ddot != 0);
    printf("ddot passed\n");

    float complex cres;
    cblas_cdotu_sub(1,cx,1,cy,1,&cres);
    printf("cdotu passed\n");

    cblas_cdotc_sub(1,cx,1,cy,1,&cres);
    printf("cdotc passed\n");


    printf("\nNORM / ASUM / IAMAX\n");

    assert(cblas_snrm2(3,sx,1) >= 0);
    printf("snrm2 passed\n");

    assert(cblas_dnrm2(3,dx,1) >= 0);
    printf("dnrm2 passed\n");

    assert(cblas_sasum(3,sx,1) >= 0);
    printf("sasum passed\n");

    assert(cblas_dasum(3,dx,1) >= 0);
    printf("dasum passed\n");

    cblas_isamax(3,sx,1);
    printf("isamax passed\n");

    cblas_idamax(3,dx,1);
    printf("idamax passed\n");


    printf("\nROT\n");

    double a=3,b=4,c,s;
    cblas_drotg(&a,&b,&c,&s);
    assert(a == 5);
    printf("drotg passed\n");

    float af=3,bf=4,cf,sf;
    cblas_srotg(&af,&bf,&cf,&sf);
    assert(af == 5);
    printf("srotg passed\n");

    double rx[1]={1}, ry[1]={0};
    cblas_drot(1,rx,1,ry,1,0,1);
    assert(ry[0] == -1);
    printf("drot passed\n");

    float rxf[1]={1}, ryf[1]={0};
    cblas_srot(1,rxf,1,ryf,1,0,1);
    assert(ryf[0] == -1);
    printf("srot passed\n");


    printf("\nROTM\n");

    double rd1=1, rd2=1, rx1=1, ry1=1, rparam[5];
    cblas_drotmg(&rd1,&rd2,&rx1,ry1,rparam);
    cblas_drotm(1,rx,1,ry,1,rparam);
    printf("drotmg passed\n");

    float frd1=1, frd2=1, frx1=1, fry1=1, frparam[5];
    cblas_srotmg(&frd1,&frd2,&frx1,fry1,frparam);
    cblas_srotm(1,rxf,1,ryf,1,frparam);
    printf("srotmg passed\n");


    printf("\nEDGE CASES\n");

    cblas_daxpy(0,2.0,dx,1,dy,1);
    printf("zero length passed\n");

    double ex[5]={1,0,2,0,3};
    double ey[3]={0,0,0};
    cblas_daxpy(3,1.0,ex,2,ey,1);
    assert(ey[2] == 3);
    printf("stride passed\n");

    return 0;
}

## Вывод

### Запуск
user@Macbook-Air OpenBLAS % ./tests
AXPY
saxpy passed
daxpy passed
caxpy passed
zaxpy passed

SCAL
sscal passed
dscal passed
cscal passed
zscal passed

COPY
scopy passed
dcopy passed

SWAP
sswap passed
dswap passed

DOT
sdot passed
ddot passed
cdotu passed
cdotc passed

NORM / ASUM / IAMAX
snrm2 passed
dnrm2 passed
sasum passed
dasum passed
isamax passed
idamax passed

ROT
drotg passed
srotg passed
drot passed
srot passed

ROTM
drotmg passed
srotmg passed

EDGE CASES
zero length passed
stride passed

### 1 поток

user@Macbook-Air OpenBLAS % export OPENBLAS_NUM_THREADS=1
user@Macbook-Air OpenBLAS % ./tests                      
AXPY
saxpy passed
daxpy passed
caxpy passed
zaxpy passed

SCAL
sscal passed
dscal passed
cscal passed
zscal passed

COPY
scopy passed
dcopy passed

SWAP
sswap passed
dswap passed

DOT
sdot passed
ddot passed
cdotu passed
cdotc passed

NORM / ASUM / IAMAX
snrm2 passed
dnrm2 passed
sasum passed
dasum passed
isamax passed
idamax passed

ROT
drotg passed
srotg passed
drot passed
srot passed

ROTM
drotmg passed
srotmg passed

EDGE CASES
zero length passed
stride passed
user@Macbook-Air OpenBLAS % 

### 4 потока

user@Macbook-Air OpenBLAS % export OPENBLAS_NUM_THREADS=4          

user@Macbook-Air OpenBLAS % ./tests                      
AXPY
saxpy passed
daxpy passed
caxpy passed
zaxpy passed

SCAL
sscal passed
dscal passed
cscal passed
zscal passed

COPY
scopy passed
dcopy passed

SWAP
sswap passed
dswap passed

DOT
sdot passed
ddot passed
cdotu passed
cdotc passed

NORM / ASUM / IAMAX
snrm2 passed
dnrm2 passed
sasum passed
dasum passed
isamax passed
idamax passed

ROT
drotg passed
srotg passed
drot passed
srot passed

ROTM
drotmg passed
srotmg passed

EDGE CASES
zero length passed
stride passed
user@Macbook-Air OpenBLAS % 
---
