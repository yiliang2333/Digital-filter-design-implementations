#ifndef __BUTTER_C
#define __BUTTER_C

/* 
    all filter coeff will be returned as: a0,a1,a2,-b1,-b2 
    b1, b2 added '-' to calculate MAC 5 times
*/

#include "butter.h"
#include "complex.h"

// PI is defined in "arm_math.h"
#ifndef PI
    #define PI               3.14159265358979f
#endif

/* 
    one-order to eight-order butterworth low-pass and high-pass filter design 
    return: == 0, OK
            != 0, ERROR,    1, order is out of range
*/
int FilterDesignButterLPF(int order, float fc, float fs, float *coeff)
{
    if ((order > 8) || (order < 1))
    {
        return 1;
    }
    fc = (fc < 20.0f) ? 20.0f : fc;
    fc = (fc > 20000.0f) ? 20000.0f : fc;
    fc = (fc > (0.49f * fs)) ? (0.49f * fs) : fc;
    
    int orderMax = 5;
    if (fc >= 50.0f)
    {
        orderMax = 6;
    }
    if (fc >= 100.0f)
    {
        orderMax = 7;
    }
    if (fc >= 200.0f)
    {
        orderMax = 8;
    }
    order = (order > orderMax) ? orderMax : order;
    int index;
    COMPLEX ct1, ct2;
    COMPLEX pl[8];
    for (index = 0; index < order; index++)
    {
        float t = (float)(index * 2 + 1) * PI / 2.0f / (float)order;
        float tt = tanf(PI * fc / fs) * fs * 2.0f;
        ct1.real = tt * (-sinf(t));
        ct1.imag = tt * (cosf(t));
        ct2.real = fs * 2.0f;
        ct2.imag = 0;
        COMPLEX pre = CDivC(ct1, ct2);
        ct1.real = 1.0f + pre.real;
        ct1.imag = pre.imag;
        ct2.real = 1.0f - pre.real;
        ct2.imag = -pre.imag;
        pl[index] = CDivC(ct1, ct2);
    }
    index = 0;
    while (index + 1 < order)
    {
        float t1 = pl[index / 2].real * pl[index / 2].real + pl[index / 2].imag * pl[index / 2].imag;
        float t2 = -2.0f * pl[index / 2].real;
        float g = (t1 + t2 + 1.0f) / 4.0f;
        *coeff++ = g;
        *coeff++ = g * 2.0f;
        *coeff++ = g;
        *coeff++ = -t2;
        *coeff++ = -t1;
        index += 2;
    }
    if (index < order)
    {
        float t2 = - pl[index / 2].real;
        float g = (1.0f + t2) / 2.0f;
        *coeff++ = g;
        *coeff++ = g;
        *coeff++ = 0;
        *coeff++ = -t2;
        *coeff++ = 0;
    }
    return 0;
}

/* 
    one-order to eight-order butterworth low-pass and high-pass filter design 
    return: == 0, OK
            != 0, ERROR,    1, order is out of range
*/
int FilterDesignButterHPF(int order, float fc, float fs, float *coeff)
{
    if ((order > 8) || (order < 1))
    {
        return 1;
    }
    fc = (fc < 20.0f) ? 20.0f : fc;
    fc = (fc > 20000.0f) ? 20000.0f : fc;
    fc = (fc > (0.49f * fs)) ? (0.49f * fs) : fc;
    int orderMax = 5;
    if (fc >= 50.0f)
    {
        orderMax = 6;
    }
    if (fc >= 100.0f)
    {
        orderMax = 7;
    }
    if (fc >= 200.0f)
    {
        orderMax = 8;
    }
    order = (order > orderMax) ? orderMax : order;
    int index;
    COMPLEX ct1, ct2;
    COMPLEX pl[8];
    for (index = 0; index < order; index++)
    {
        float t = (float)(index * 2 + 1) * PI / 2.0f / (float)order;
        float tt = tanf(PI * fc / fs) * fs * 2.0f;
        ct1.real = tt;
        ct1.imag = 0;
        ct2.real = -sinf(t);
        ct2.imag =  cosf(t);
        ct1 = CDivC(ct1, ct2);
        ct2.real = fs * 2.0f;
        ct2.imag = 0;
        COMPLEX pre = CDivC(ct1, ct2);
        ct1.real = 1.0f + pre.real;
        ct1.imag = pre.imag;
        ct2.real = 1.0f - pre.real;
        ct2.imag = -pre.imag;
        pl[index] = CDivC(ct1, ct2);
    }
    index = 0;
    while (index + 1 < order)
    {
        float t1 = pl[index / 2].real * pl[index / 2].real + pl[index / 2].imag * pl[index / 2].imag;
        float t2 = -2.0f * pl[index / 2].real;
        float g = (t1 + fabsf(t2) + 1.0f) / 4.0f;
        *coeff++ = g;
        *coeff++ = g * -2.0f;
        *coeff++ = g;
        *coeff++ = -t2;
        *coeff++ = -t1;
        index += 2;
    }
    if (index < order)
    {
        float t2 = - pl[index / 2].real;
        float g = (1.0f + fabsf(t2)) / 2.0f;
        *coeff++ = g;
        *coeff++ = g * -1.0f;
        *coeff++ = 0;
        *coeff++ = -t2;
        *coeff++ = 0;
    }
    return 0;
}

/* 
    one-order and two-order butterworth low-pass and high-pass filter design (type == 0, low-pass, type == 1, high-pass)
    return: == 0, OK
            != 0, ERROR,    1, order is out of range
                            2, type is out of range
*/
int FilterDesignButter(int order, float fc, float fs, float *coeff, int type)
{
    if ((order != 1) || (order != 2))
    {
        return 1;
    }
    if ((type != 0) || (type != 1))
    {
        return 2;
    }
    fc = (fc < 5.0f) ? 5.0f : fc;
    float wn = fc / fs;
    wn = (wn > 0.49f) ? 0.49f : wn;
    float a0, a1, a2, b1, b2;
    a0 = 1.0f;
    a1 = a2 = b1 = b2 = 0;
    if (order == 1)
    {
        float cTanf = tanf(PI * wn);
        float c = (cTanf - 1.0f) / (cTanf + 1.0f);
        b1 = c;
        if (type == 0)
        {
            a0 = (1.0f + c) / 2.0f;
            a1 = a0;
        }
        else //if (type == 1)   // judge the value of type at the beginnig, type only can be '0' or '1', so here do not need to judge again
        {
            a0 = (1.0f - c) / 2.0f;
            a1 = (c - 1.0f) / 2.0f;
        }
    }
    else //if (order == 2)      // judge the value of order at the beginnig, type only can be '1' or '2', so here do not need to judge again
    {
        float K = tanf(PI * wn);
        float t = sqrtf(2.0f) * K + K * K;
        float c = 1.0f + t;
        b1 = 2.0f * (K * K - 1.0f) / c;
        b2 = (1.0f - t) / c;
        if (type == 0)
        {
            a0 = K * K / c;
            a1 = 2.0f * a0;
            a2 = a0;
        }
        else //if (type == 1)   // judge the value of type at the beginnig, type only can be '0' or '1', so here do not need to judge again
        {
            a0 = 1.0f / c;
            a1 = -2.0f * a0;
            a2 = a0;
        }
    }
    *coeff++ = a0;
    *coeff++ = a1;
    *coeff++ = a2;
    *coeff++ = b1;
    *coeff++ = b2;
    return 0;
}

#endif




