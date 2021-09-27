  /* base */
#define HZ                16
#define STEPS             128
#define BUFFER            16
#define TUNE              0
#define CHANNELS          2
#define FFT               0 // 0 -> Raw Data, 1 -> 1D FFT, 2 -> 2D FFT
  /* 1d fft plot */
#define FFT1_PLOT_TYPE    0 // 0 -> 2D Plot, 1 -> 3D Plot
#define FFT1_MAG_MAX      2
#define FFT1_MAG_MIN      0
  /* 2d fft plot */
#define FFT2_PLOT_TYPE    0 // 0 -> 2D Plot, 1 -> 3D Plot
#define FFT2_MAG_MAX      2
#define FFT2_MAG_MIN      0
#define FFT2_SUM_MAX      1000
#define FFT2_SUM_MIN      0
#define FFT2_CT_MAX       1000
#define CROSSTALK_SAMPLES 1
  /* data processing */
#define ZERO_PAD          512
#define DC_OFFSET         true
#define HAMMING_WINDOW    true
#define ZERO_FREE_SPACE   false /* CAUTION! */
  /* callibration */
#define FREE_SPACE_CYCLES 0 /* CAUTION! */
#define RADAR_CYCLES      0 /* CAUTION! */
#define FFT_CYCLES        0 /* CAUTION! */
