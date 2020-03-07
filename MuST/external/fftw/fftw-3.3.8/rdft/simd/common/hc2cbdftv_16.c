/*
 * Copyright (c) 2003, 2007-14 Matteo Frigo
 * Copyright (c) 2003, 2007-14 Massachusetts Institute of Technology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

/* This file was automatically generated --- DO NOT EDIT */
/* Generated on Thu May 24 08:08:12 EDT 2018 */

#include "rdft/codelet-rdft.h"

#if defined(ARCH_PREFERS_FMA) || defined(ISA_EXTENSION_PREFERS_FMA)

/* Generated by: ../../../genfft/gen_hc2cdft_c.native -fma -simd -compact -variables 4 -pipeline-latency 8 -trivial-stores -variables 32 -no-generate-bytw -n 16 -dif -sign 1 -name hc2cbdftv_16 -include rdft/simd/hc2cbv.h */

/*
 * This function contains 103 FP additions, 80 FP multiplications,
 * (or, 53 additions, 30 multiplications, 50 fused multiply/add),
 * 79 stack variables, 3 constants, and 32 memory accesses
 */
#include "rdft/simd/hc2cbv.h"

static void hc2cbdftv_16(R *Rp, R *Ip, R *Rm, R *Im, const R *W, stride rs, INT mb, INT me, INT ms)
{
     DVK(KP923879532, +0.923879532511286756128183189396788286822416626);
     DVK(KP414213562, +0.414213562373095048801688724209698078569671875);
     DVK(KP707106781, +0.707106781186547524400844362104849039284835938);
     {
	  INT m;
	  for (m = mb, W = W + ((mb - 1) * ((TWVL / VL) * 30)); m < me; m = m + VL, Rp = Rp + (VL * ms), Ip = Ip + (VL * ms), Rm = Rm - (VL * ms), Im = Im - (VL * ms), W = W + (TWVL * 30), MAKE_VOLATILE_STRIDE(64, rs)) {
	       V T8, Tv, TE, T1t, TP, T1w, T10, T1p, Tn, Tw, T13, T1q, TL, T1x, TS;
	       V T1u;
	       {
		    V T4, TA, Tu, TC, T7, TN, Tr, TB, T2, T3, Ts, Tt, T5, T6, Tp;
		    V Tq, TD, TO, TY, TZ, Tb, TF, Tl, TJ, Te, TG, Ti, TI, T9, Ta;
		    V Tj, Tk, Tc, Td, Tg, Th, Tf, Tm, T11, T12, TH, TK, TQ, TR;
		    T2 = LD(&(Rp[0]), ms, &(Rp[0]));
		    T3 = LD(&(Rm[WS(rs, 7)]), -ms, &(Rm[WS(rs, 1)]));
		    T4 = VFMACONJ(T3, T2);
		    TA = VFNMSCONJ(T3, T2);
		    Ts = LD(&(Rp[WS(rs, 6)]), ms, &(Rp[0]));
		    Tt = LD(&(Rm[WS(rs, 1)]), -ms, &(Rm[WS(rs, 1)]));
		    Tu = VFMACONJ(Tt, Ts);
		    TC = VFMSCONJ(Tt, Ts);
		    T5 = LD(&(Rp[WS(rs, 4)]), ms, &(Rp[0]));
		    T6 = LD(&(Rm[WS(rs, 3)]), -ms, &(Rm[WS(rs, 1)]));
		    T7 = VFMACONJ(T6, T5);
		    TN = VFNMSCONJ(T6, T5);
		    Tp = LD(&(Rp[WS(rs, 2)]), ms, &(Rp[0]));
		    Tq = LD(&(Rm[WS(rs, 5)]), -ms, &(Rm[WS(rs, 1)]));
		    Tr = VFMACONJ(Tq, Tp);
		    TB = VFNMSCONJ(Tq, Tp);
		    T8 = VSUB(T4, T7);
		    Tv = VSUB(Tr, Tu);
		    TD = VADD(TB, TC);
		    TE = VFMA(LDK(KP707106781), TD, TA);
		    T1t = VFNMS(LDK(KP707106781), TD, TA);
		    TO = VSUB(TB, TC);
		    TP = VFMA(LDK(KP707106781), TO, TN);
		    T1w = VFNMS(LDK(KP707106781), TO, TN);
		    TY = VADD(T4, T7);
		    TZ = VADD(Tr, Tu);
		    T10 = VADD(TY, TZ);
		    T1p = VSUB(TY, TZ);
		    T9 = LD(&(Rp[WS(rs, 1)]), ms, &(Rp[WS(rs, 1)]));
		    Ta = LD(&(Rm[WS(rs, 6)]), -ms, &(Rm[0]));
		    Tb = VFMACONJ(Ta, T9);
		    TF = VFNMSCONJ(Ta, T9);
		    Tj = LD(&(Rp[WS(rs, 3)]), ms, &(Rp[WS(rs, 1)]));
		    Tk = LD(&(Rm[WS(rs, 4)]), -ms, &(Rm[0]));
		    Tl = VFMACONJ(Tk, Tj);
		    TJ = VFNMSCONJ(Tk, Tj);
		    Tc = LD(&(Rp[WS(rs, 5)]), ms, &(Rp[WS(rs, 1)]));
		    Td = LD(&(Rm[WS(rs, 2)]), -ms, &(Rm[0]));
		    Te = VFMACONJ(Td, Tc);
		    TG = VFNMSCONJ(Td, Tc);
		    Tg = LD(&(Rp[WS(rs, 7)]), ms, &(Rp[WS(rs, 1)]));
		    Th = LD(&(Rm[0]), -ms, &(Rm[0]));
		    Ti = VFMACONJ(Th, Tg);
		    TI = VFMSCONJ(Th, Tg);
		    Tf = VSUB(Tb, Te);
		    Tm = VSUB(Ti, Tl);
		    Tn = VADD(Tf, Tm);
		    Tw = VSUB(Tf, Tm);
		    T11 = VADD(Tb, Te);
		    T12 = VADD(Ti, Tl);
		    T13 = VADD(T11, T12);
		    T1q = VSUB(T11, T12);
		    TH = VFNMS(LDK(KP414213562), TG, TF);
		    TK = VFMA(LDK(KP414213562), TJ, TI);
		    TL = VADD(TH, TK);
		    T1x = VSUB(TH, TK);
		    TQ = VFMA(LDK(KP414213562), TF, TG);
		    TR = VFNMS(LDK(KP414213562), TI, TJ);
		    TS = VADD(TQ, TR);
		    T1u = VSUB(TQ, TR);
	       }
	       {
		    V T1j, T1R, T1c, T1J, T1g, T1l, T1N, T1T, T1Q, T1a, T1b, T19, T1I, T1e, T1f;
		    V T1d, T1k, T1L, T1M, T1K, T1S, T1h, T1U, T1V, T1i, T1m, T1O, T1P, T1n, T14;
		    V T1r, Ty, T1D, TU, T16, T1z, T1F, TX, T1o, To, Tx, T1, T1C, TM, TT;
		    V Tz, T15, T1v, T1y, T1s, T1E, TV, T1G, T1H, TW, T17, T1A, T1B, T18;
		    T1j = VADD(T10, T13);
		    T1Q = LDW(&(W[TWVL * 22]));
		    T1R = VZMUL(T1Q, VFNMSI(T1q, T1p));
		    T1a = VFMA(LDK(KP707106781), Tn, T8);
		    T1b = VFMA(LDK(KP707106781), Tw, Tv);
		    T19 = LDW(&(W[TWVL * 26]));
		    T1c = VZMUL(T19, VFNMSI(T1b, T1a));
		    T1I = LDW(&(W[TWVL * 2]));
		    T1J = VZMUL(T1I, VFMAI(T1b, T1a));
		    T1e = VFMA(LDK(KP923879532), TL, TE);
		    T1f = VFMA(LDK(KP923879532), TS, TP);
		    T1d = LDW(&(W[TWVL * 28]));
		    T1g = VZMULI(T1d, VFNMSI(T1f, T1e));
		    T1k = LDW(&(W[0]));
		    T1l = VZMULI(T1k, VFMAI(T1f, T1e));
		    T1L = VFMA(LDK(KP923879532), T1u, T1t);
		    T1M = VFNMS(LDK(KP923879532), T1x, T1w);
		    T1K = LDW(&(W[TWVL * 4]));
		    T1N = VZMULI(T1K, VFNMSI(T1M, T1L));
		    T1S = LDW(&(W[TWVL * 24]));
		    T1T = VZMULI(T1S, VFMAI(T1M, T1L));
		    T1h = VCONJ(VSUB(T1c, T1g));
		    ST(&(Rm[WS(rs, 7)]), T1h, -ms, &(Rm[WS(rs, 1)]));
		    T1U = VCONJ(VSUB(T1R, T1T));
		    ST(&(Rm[WS(rs, 6)]), T1U, -ms, &(Rm[0]));
		    T1V = VADD(T1R, T1T);
		    ST(&(Rp[WS(rs, 6)]), T1V, ms, &(Rp[0]));
		    T1i = VADD(T1c, T1g);
		    ST(&(Rp[WS(rs, 7)]), T1i, ms, &(Rp[WS(rs, 1)]));
		    T1m = VCONJ(VSUB(T1j, T1l));
		    ST(&(Rm[0]), T1m, -ms, &(Rm[0]));
		    T1O = VCONJ(VSUB(T1J, T1N));
		    ST(&(Rm[WS(rs, 1)]), T1O, -ms, &(Rm[WS(rs, 1)]));
		    T1P = VADD(T1J, T1N);
		    ST(&(Rp[WS(rs, 1)]), T1P, ms, &(Rp[WS(rs, 1)]));
		    T1n = VADD(T1j, T1l);
		    ST(&(Rp[0]), T1n, ms, &(Rp[0]));
		    TX = LDW(&(W[TWVL * 14]));
		    T14 = VZMUL(TX, VSUB(T10, T13));
		    T1o = LDW(&(W[TWVL * 6]));
		    T1r = VZMUL(T1o, VFMAI(T1q, T1p));
		    To = VFNMS(LDK(KP707106781), Tn, T8);
		    Tx = VFNMS(LDK(KP707106781), Tw, Tv);
		    T1 = LDW(&(W[TWVL * 10]));
		    Ty = VZMUL(T1, VFNMSI(Tx, To));
		    T1C = LDW(&(W[TWVL * 18]));
		    T1D = VZMUL(T1C, VFMAI(Tx, To));
		    TM = VFNMS(LDK(KP923879532), TL, TE);
		    TT = VFNMS(LDK(KP923879532), TS, TP);
		    Tz = LDW(&(W[TWVL * 12]));
		    TU = VZMULI(Tz, VFNMSI(TT, TM));
		    T15 = LDW(&(W[TWVL * 16]));
		    T16 = VZMULI(T15, VFMAI(TT, TM));
		    T1v = VFNMS(LDK(KP923879532), T1u, T1t);
		    T1y = VFMA(LDK(KP923879532), T1x, T1w);
		    T1s = LDW(&(W[TWVL * 8]));
		    T1z = VZMULI(T1s, VFMAI(T1y, T1v));
		    T1E = LDW(&(W[TWVL * 20]));
		    T1F = VZMULI(T1E, VFNMSI(T1y, T1v));
		    TV = VCONJ(VSUB(Ty, TU));
		    ST(&(Rm[WS(rs, 3)]), TV, -ms, &(Rm[WS(rs, 1)]));
		    T1G = VCONJ(VSUB(T1D, T1F));
		    ST(&(Rm[WS(rs, 5)]), T1G, -ms, &(Rm[WS(rs, 1)]));
		    T1H = VADD(T1D, T1F);
		    ST(&(Rp[WS(rs, 5)]), T1H, ms, &(Rp[WS(rs, 1)]));
		    TW = VADD(Ty, TU);
		    ST(&(Rp[WS(rs, 3)]), TW, ms, &(Rp[WS(rs, 1)]));
		    T17 = VCONJ(VSUB(T14, T16));
		    ST(&(Rm[WS(rs, 4)]), T17, -ms, &(Rm[0]));
		    T1A = VCONJ(VSUB(T1r, T1z));
		    ST(&(Rm[WS(rs, 2)]), T1A, -ms, &(Rm[0]));
		    T1B = VADD(T1r, T1z);
		    ST(&(Rp[WS(rs, 2)]), T1B, ms, &(Rp[0]));
		    T18 = VADD(T14, T16);
		    ST(&(Rp[WS(rs, 4)]), T18, ms, &(Rp[0]));
	       }
	  }
     }
     VLEAVE();
}

static const tw_instr twinstr[] = {
     VTW(1, 1),
     VTW(1, 2),
     VTW(1, 3),
     VTW(1, 4),
     VTW(1, 5),
     VTW(1, 6),
     VTW(1, 7),
     VTW(1, 8),
     VTW(1, 9),
     VTW(1, 10),
     VTW(1, 11),
     VTW(1, 12),
     VTW(1, 13),
     VTW(1, 14),
     VTW(1, 15),
     {TW_NEXT, VL, 0}
};

static const hc2c_desc desc = { 16, XSIMD_STRING("hc2cbdftv_16"), twinstr, &GENUS, {53, 30, 50, 0} };

void XSIMD(codelet_hc2cbdftv_16) (planner *p) {
     X(khc2c_register) (p, hc2cbdftv_16, &desc, HC2C_VIA_DFT);
}
#else

/* Generated by: ../../../genfft/gen_hc2cdft_c.native -simd -compact -variables 4 -pipeline-latency 8 -trivial-stores -variables 32 -no-generate-bytw -n 16 -dif -sign 1 -name hc2cbdftv_16 -include rdft/simd/hc2cbv.h */

/*
 * This function contains 103 FP additions, 42 FP multiplications,
 * (or, 99 additions, 38 multiplications, 4 fused multiply/add),
 * 83 stack variables, 3 constants, and 32 memory accesses
 */
#include "rdft/simd/hc2cbv.h"

static void hc2cbdftv_16(R *Rp, R *Ip, R *Rm, R *Im, const R *W, stride rs, INT mb, INT me, INT ms)
{
     DVK(KP382683432, +0.382683432365089771728459984030398866761344562);
     DVK(KP923879532, +0.923879532511286756128183189396788286822416626);
     DVK(KP707106781, +0.707106781186547524400844362104849039284835938);
     {
	  INT m;
	  for (m = mb, W = W + ((mb - 1) * ((TWVL / VL) * 30)); m < me; m = m + VL, Rp = Rp + (VL * ms), Ip = Ip + (VL * ms), Rm = Rm - (VL * ms), Im = Im - (VL * ms), W = W + (TWVL * 30), MAKE_VOLATILE_STRIDE(64, rs)) {
	       V Tf, T16, TZ, T1C, TI, T1a, TV, T1D, T1F, T1G, Ty, T19, TC, T17, TS;
	       V T10;
	       {
		    V T2, TD, T4, TF, Tc, Tb, Td, T6, T8, T9, T3, TE, Ta, T7, T5;
		    V Te, TX, TY, TG, TH, TT, TU, Tj, TM, Tw, TQ, Tn, TN, Ts, TP;
		    V Tg, Ti, Th, Tt, Tv, Tu, Tk, Tm, Tl, Tr, Tq, Tp, To, Tx, TA;
		    V TB, TO, TR;
		    T2 = LD(&(Rp[0]), ms, &(Rp[0]));
		    TD = LD(&(Rp[WS(rs, 4)]), ms, &(Rp[0]));
		    T3 = LD(&(Rm[WS(rs, 7)]), -ms, &(Rm[WS(rs, 1)]));
		    T4 = VCONJ(T3);
		    TE = LD(&(Rm[WS(rs, 3)]), -ms, &(Rm[WS(rs, 1)]));
		    TF = VCONJ(TE);
		    Tc = LD(&(Rp[WS(rs, 6)]), ms, &(Rp[0]));
		    Ta = LD(&(Rm[WS(rs, 1)]), -ms, &(Rm[WS(rs, 1)]));
		    Tb = VCONJ(Ta);
		    Td = VSUB(Tb, Tc);
		    T6 = LD(&(Rp[WS(rs, 2)]), ms, &(Rp[0]));
		    T7 = LD(&(Rm[WS(rs, 5)]), -ms, &(Rm[WS(rs, 1)]));
		    T8 = VCONJ(T7);
		    T9 = VSUB(T6, T8);
		    T5 = VSUB(T2, T4);
		    Te = VMUL(LDK(KP707106781), VADD(T9, Td));
		    Tf = VADD(T5, Te);
		    T16 = VSUB(T5, Te);
		    TX = VADD(T2, T4);
		    TY = VADD(TD, TF);
		    TZ = VSUB(TX, TY);
		    T1C = VADD(TX, TY);
		    TG = VSUB(TD, TF);
		    TH = VMUL(LDK(KP707106781), VSUB(T9, Td));
		    TI = VADD(TG, TH);
		    T1a = VSUB(TH, TG);
		    TT = VADD(T6, T8);
		    TU = VADD(Tb, Tc);
		    TV = VSUB(TT, TU);
		    T1D = VADD(TT, TU);
		    Tg = LD(&(Rp[WS(rs, 1)]), ms, &(Rp[WS(rs, 1)]));
		    Th = LD(&(Rm[WS(rs, 6)]), -ms, &(Rm[0]));
		    Ti = VCONJ(Th);
		    Tj = VSUB(Tg, Ti);
		    TM = VADD(Tg, Ti);
		    Tt = LD(&(Rp[WS(rs, 3)]), ms, &(Rp[WS(rs, 1)]));
		    Tu = LD(&(Rm[WS(rs, 4)]), -ms, &(Rm[0]));
		    Tv = VCONJ(Tu);
		    Tw = VSUB(Tt, Tv);
		    TQ = VADD(Tt, Tv);
		    Tk = LD(&(Rp[WS(rs, 5)]), ms, &(Rp[WS(rs, 1)]));
		    Tl = LD(&(Rm[WS(rs, 2)]), -ms, &(Rm[0]));
		    Tm = VCONJ(Tl);
		    Tn = VSUB(Tk, Tm);
		    TN = VADD(Tk, Tm);
		    Tr = LD(&(Rp[WS(rs, 7)]), ms, &(Rp[WS(rs, 1)]));
		    Tp = LD(&(Rm[0]), -ms, &(Rm[0]));
		    Tq = VCONJ(Tp);
		    Ts = VSUB(Tq, Tr);
		    TP = VADD(Tq, Tr);
		    T1F = VADD(TM, TN);
		    T1G = VADD(TP, TQ);
		    To = VFNMS(LDK(KP382683432), Tn, VMUL(LDK(KP923879532), Tj));
		    Tx = VFMA(LDK(KP923879532), Ts, VMUL(LDK(KP382683432), Tw));
		    Ty = VADD(To, Tx);
		    T19 = VSUB(To, Tx);
		    TA = VFMA(LDK(KP382683432), Tj, VMUL(LDK(KP923879532), Tn));
		    TB = VFNMS(LDK(KP382683432), Ts, VMUL(LDK(KP923879532), Tw));
		    TC = VADD(TA, TB);
		    T17 = VSUB(TA, TB);
		    TO = VSUB(TM, TN);
		    TR = VSUB(TP, TQ);
		    TS = VMUL(LDK(KP707106781), VSUB(TO, TR));
		    T10 = VMUL(LDK(KP707106781), VADD(TO, TR));
	       }
	       {
		    V T21, T1W, T1u, T20, T1I, T1O, TK, T1S, T12, T1e, T1k, T1A, T1o, T1w, T1c;
		    V T1M, T1U, T1V, T1T, T1s, T1t, T1r, T1Z, T1E, T1H, T1B, T1N, Tz, TJ, T1;
		    V T1R, TW, T11, TL, T1d, T1i, T1j, T1h, T1z, T1m, T1n, T1l, T1v, T18, T1b;
		    V T15, T1L, T13, T1g, T1X, T23, T14, T1f, T1Y, T22, T1p, T1y, T1J, T1Q, T1q;
		    V T1x, T1K, T1P;
		    T1U = VADD(T1C, T1D);
		    T1V = VADD(T1F, T1G);
		    T21 = VADD(T1U, T1V);
		    T1T = LDW(&(W[TWVL * 14]));
		    T1W = VZMUL(T1T, VSUB(T1U, T1V));
		    T1s = VADD(Tf, Ty);
		    T1t = VBYI(VADD(TI, TC));
		    T1r = LDW(&(W[TWVL * 28]));
		    T1u = VZMULI(T1r, VSUB(T1s, T1t));
		    T1Z = LDW(&(W[0]));
		    T20 = VZMULI(T1Z, VADD(T1s, T1t));
		    T1E = VSUB(T1C, T1D);
		    T1H = VBYI(VSUB(T1F, T1G));
		    T1B = LDW(&(W[TWVL * 22]));
		    T1I = VZMUL(T1B, VSUB(T1E, T1H));
		    T1N = LDW(&(W[TWVL * 6]));
		    T1O = VZMUL(T1N, VADD(T1E, T1H));
		    Tz = VSUB(Tf, Ty);
		    TJ = VBYI(VSUB(TC, TI));
		    T1 = LDW(&(W[TWVL * 12]));
		    TK = VZMULI(T1, VADD(Tz, TJ));
		    T1R = LDW(&(W[TWVL * 16]));
		    T1S = VZMULI(T1R, VSUB(Tz, TJ));
		    TW = VBYI(VSUB(TS, TV));
		    T11 = VSUB(TZ, T10);
		    TL = LDW(&(W[TWVL * 10]));
		    T12 = VZMUL(TL, VADD(TW, T11));
		    T1d = LDW(&(W[TWVL * 18]));
		    T1e = VZMUL(T1d, VSUB(T11, TW));
		    T1i = VBYI(VADD(T1a, T19));
		    T1j = VADD(T16, T17);
		    T1h = LDW(&(W[TWVL * 4]));
		    T1k = VZMULI(T1h, VADD(T1i, T1j));
		    T1z = LDW(&(W[TWVL * 24]));
		    T1A = VZMULI(T1z, VSUB(T1j, T1i));
		    T1m = VBYI(VADD(TV, TS));
		    T1n = VADD(TZ, T10);
		    T1l = LDW(&(W[TWVL * 2]));
		    T1o = VZMUL(T1l, VADD(T1m, T1n));
		    T1v = LDW(&(W[TWVL * 26]));
		    T1w = VZMUL(T1v, VSUB(T1n, T1m));
		    T18 = VSUB(T16, T17);
		    T1b = VBYI(VSUB(T19, T1a));
		    T15 = LDW(&(W[TWVL * 20]));
		    T1c = VZMULI(T15, VSUB(T18, T1b));
		    T1L = LDW(&(W[TWVL * 8]));
		    T1M = VZMULI(T1L, VADD(T1b, T18));
		    T13 = VADD(TK, T12);
		    ST(&(Rp[WS(rs, 3)]), T13, ms, &(Rp[WS(rs, 1)]));
		    T1g = VCONJ(VSUB(T1e, T1c));
		    ST(&(Rm[WS(rs, 5)]), T1g, -ms, &(Rm[WS(rs, 1)]));
		    T1X = VADD(T1S, T1W);
		    ST(&(Rp[WS(rs, 4)]), T1X, ms, &(Rp[0]));
		    T23 = VCONJ(VSUB(T21, T20));
		    ST(&(Rm[0]), T23, -ms, &(Rm[0]));
		    T14 = VCONJ(VSUB(T12, TK));
		    ST(&(Rm[WS(rs, 3)]), T14, -ms, &(Rm[WS(rs, 1)]));
		    T1f = VADD(T1c, T1e);
		    ST(&(Rp[WS(rs, 5)]), T1f, ms, &(Rp[WS(rs, 1)]));
		    T1Y = VCONJ(VSUB(T1W, T1S));
		    ST(&(Rm[WS(rs, 4)]), T1Y, -ms, &(Rm[0]));
		    T22 = VADD(T20, T21);
		    ST(&(Rp[0]), T22, ms, &(Rp[0]));
		    T1p = VADD(T1k, T1o);
		    ST(&(Rp[WS(rs, 1)]), T1p, ms, &(Rp[WS(rs, 1)]));
		    T1y = VCONJ(VSUB(T1w, T1u));
		    ST(&(Rm[WS(rs, 7)]), T1y, -ms, &(Rm[WS(rs, 1)]));
		    T1J = VADD(T1A, T1I);
		    ST(&(Rp[WS(rs, 6)]), T1J, ms, &(Rp[0]));
		    T1Q = VCONJ(VSUB(T1O, T1M));
		    ST(&(Rm[WS(rs, 2)]), T1Q, -ms, &(Rm[0]));
		    T1q = VCONJ(VSUB(T1o, T1k));
		    ST(&(Rm[WS(rs, 1)]), T1q, -ms, &(Rm[WS(rs, 1)]));
		    T1x = VADD(T1u, T1w);
		    ST(&(Rp[WS(rs, 7)]), T1x, ms, &(Rp[WS(rs, 1)]));
		    T1K = VCONJ(VSUB(T1I, T1A));
		    ST(&(Rm[WS(rs, 6)]), T1K, -ms, &(Rm[0]));
		    T1P = VADD(T1M, T1O);
		    ST(&(Rp[WS(rs, 2)]), T1P, ms, &(Rp[0]));
	       }
	  }
     }
     VLEAVE();
}

static const tw_instr twinstr[] = {
     VTW(1, 1),
     VTW(1, 2),
     VTW(1, 3),
     VTW(1, 4),
     VTW(1, 5),
     VTW(1, 6),
     VTW(1, 7),
     VTW(1, 8),
     VTW(1, 9),
     VTW(1, 10),
     VTW(1, 11),
     VTW(1, 12),
     VTW(1, 13),
     VTW(1, 14),
     VTW(1, 15),
     {TW_NEXT, VL, 0}
};

static const hc2c_desc desc = { 16, XSIMD_STRING("hc2cbdftv_16"), twinstr, &GENUS, {99, 38, 4, 0} };

void XSIMD(codelet_hc2cbdftv_16) (planner *p) {
     X(khc2c_register) (p, hc2cbdftv_16, &desc, HC2C_VIA_DFT);
}
#endif
