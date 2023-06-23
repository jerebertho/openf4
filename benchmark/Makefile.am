# FLAGS

AM_CPPFLAGS= $(OPTIMISATION_FLAGS) $(DEBUG_FLAGS) $(VECTORISATION_FLAGS) -I$(top_srcdir) 

AM_LDFLAGS = 

if OPENF4_HAVE_OMP
AM_LDFLAGS+= $(OMPFLAGS)
AM_CPPFLAGS+= $(OMPFLAGS)
endif

if OPENF4_HAVE_GIVARO
AM_LDFLAGS+=$(GIVARO_LIBS) $(GMP_LIBS)
AM_CPPFLAGS+=$(GIVARO_CFLAGS) $(GMP_CFLAGS)
endif

if OPENF4_HAVE_FFLAS_FFPACK
AM_LDFLAGS+=$(CBLAS_LIBS) $(FFLAS_FFPACK_LIBS)
AM_CPPFLAGS+=$(BLAS_CFLAGS) $(FFLAS_FFPACK_CFLAGS)
endif


# Executable

noinst_PROGRAMS = benchmark_short benchmark_int benchmark_long \
benchmark_semaev benchmark_gf2 benchmark_gf2_extension_16_bits \
benchmark_gf2_extension_32_bits benchmark_gf2_extension_64_bits \
benchmark_kat9 benchmark_kat10 benchmark_kat11 benchmark_kat12 \
benchmark_kat13 benchmark_kat14 benchmark_kat15 \
benchmark_eco10 benchmark_eco11 benchmark_eco12 \
benchmark_eco13 benchmark_eco14 benchmark_eco15 \
benchmark_cp_d3_n5_p2 benchmark_cp_d3_n6_p2 \
benchmark_cp_d3_n7_p2 benchmark_cp_d3_n8_p2 \
benchmark_cp_d4_n4_p3 benchmark_cp_d4_n5_p3 \
benchmark_cp_d3_n6_p6 benchmark_cp_d4_n6_p6 benchmark_cp_d3_n7_p7 \
benchmark_sdk benchmark_all benchmark_vor benchmark_sot \
benchmark_klyz12 benchmark_klyz22 benchmark_klyz13 benchmark_klyz23 \
benchmark_noon6 benchmark_noon7 benchmark_noon8 \
benchmark_gametwo7


if OPENF4_HAVE_GIVARO
noinst_PROGRAMS += benchmark_givaro_modular_integer benchmark_givaro_gfq benchmark_givaro_modular_log16
endif


# Sources

benchmark_short_SOURCES = benchmark-short.cpp
benchmark_int_SOURCES = benchmark-int.cpp
benchmark_long_SOURCES = benchmark-long.cpp
benchmark_semaev_SOURCES = benchmark-semaev.cpp
benchmark_gf2_SOURCES = benchmark-gf2.cpp
benchmark_gf2_extension_16_bits_SOURCES = benchmark-gf2-extension-16bits.cpp
benchmark_gf2_extension_32_bits_SOURCES = benchmark-gf2-extension-32bits.cpp
benchmark_gf2_extension_64_bits_SOURCES = benchmark-gf2-extension-64bits.cpp
benchmark_kat9_SOURCES = benchmark-kat9.cpp
benchmark_kat10_SOURCES = benchmark-kat10.cpp
benchmark_kat11_SOURCES = benchmark-kat11.cpp
benchmark_kat12_SOURCES = benchmark-kat12.cpp
benchmark_kat13_SOURCES = benchmark-kat13.cpp
benchmark_kat14_SOURCES = benchmark-kat14.cpp
benchmark_kat15_SOURCES = benchmark-kat15.cpp
benchmark_eco10_SOURCES = benchmark-eco10.cpp
benchmark_eco11_SOURCES = benchmark-eco11.cpp
benchmark_eco12_SOURCES = benchmark-eco12.cpp
benchmark_eco13_SOURCES = benchmark-eco13.cpp
benchmark_eco14_SOURCES = benchmark-eco14.cpp
benchmark_eco15_SOURCES = benchmark-eco15.cpp
benchmark_cp_d3_n5_p2_SOURCES = benchmark-cp_d3_n5_p2.cpp
benchmark_cp_d3_n6_p2_SOURCES = benchmark-cp_d3_n6_p2.cpp
benchmark_cp_d3_n7_p2_SOURCES = benchmark-cp_d3_n7_p2.cpp
benchmark_cp_d3_n8_p2_SOURCES = benchmark-cp_d3_n8_p2.cpp
benchmark_cp_d4_n4_p3_SOURCES = benchmark-cp_d4_n4_p3.cpp
benchmark_cp_d4_n5_p3_SOURCES = benchmark-cp_d4_n5_p3.cpp
benchmark_cp_d3_n6_p6_SOURCES = benchmark-cp_d3_n6_p6.cpp
benchmark_cp_d4_n6_p6_SOURCES = benchmark-cp_d4_n6_p6.cpp
benchmark_cp_d3_n7_p7_SOURCES = benchmark-cp_d3_n7_p7.cpp
benchmark_sdk_SOURCES = benchmark-sdk.cpp
benchmark_all_SOURCES = benchmark-all.cpp
benchmark_vor_SOURCES = benchmark-vor.cpp
benchmark_sot_SOURCES = benchmark-sot.cpp
benchmark_klyz12_SOURCES = benchmark-klyz12.cpp
benchmark_klyz13_SOURCES = benchmark-klyz13.cpp
benchmark_klyz22_SOURCES = benchmark-klyz22.cpp
benchmark_klyz23_SOURCES = benchmark-klyz23.cpp
benchmark_noon6_SOURCES = benchmark-noon6.cpp
benchmark_noon7_SOURCES = benchmark-noon7.cpp
benchmark_noon8_SOURCES = benchmark-noon8.cpp
benchmark_gametwo7_SOURCES = benchmark-gametwo7.cpp


if OPENF4_HAVE_GIVARO
benchmark_givaro_modular_integer_SOURCES = benchmark-givaro-modular-integer.cpp
benchmark_givaro_gfq_SOURCES = benchmark-givaro-gfq.cpp
benchmark_givaro_modular_log16_SOURCES = benchmark-givaro-modular-log16.cpp
endif

