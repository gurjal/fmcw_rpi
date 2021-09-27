#include <fstream>
#include <iostream>
#include <vector>
#include <thread>
#include <math.h>
#include <bcm2835.h>
#include <fftw3.h>
#include <unistd.h>
#include <sched.h>
#include "config.h"

#define REAL 0
#define IMAG 1

#if CHANNELS==1
#define SPI_LEN 3
static char _softspan_data[SPI_LEN] = { 0x00, 0x80, 0x00 };
static char _send_data[SPI_LEN] = { 0x00, 0x00, 0x00 };
static char _read_data[SPI_LEN] = { 0x00, 0x00, 0x00 };
#elif CHANNELS==2
#define SPI_LEN 6
static char _softspan_data[SPI_LEN] = { 0x00, 0x90, 0x00, 0x00, 0x00, 0x00 };
static char _send_data[SPI_LEN] = { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };
static char _read_data[SPI_LEN] = { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };
#elif CHANNELS==3
#define SPI_LEN 9
static char _softspan_data[SPI_LEN] = { 0x04, 0x90, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };
static char _send_data[SPI_LEN] = { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };
static char _read_data[SPI_LEN] = { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };
#elif CHANNELS==4
#define SPI_LEN 12
static char _softspan_data[SPI_LEN] = { 0x24, 0x90, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };
static char _send_data[SPI_LEN] = { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };
static char _read_data[SPI_LEN] = { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };
#endif

static double dc_offset[BUFFER];
static double hamming_window[STEPS];
static long double radar_calibration[STEPS] = { 0 };
static long double free_space_calibration[STEPS] = { 0 };
static double fft_calibration[(ZERO_PAD+STEPS)/2] = { 0 };
static int ZERO_PAD_STEPS = ZERO_PAD+STEPS;
static fftw_complex* _raw_data[CHANNELS][BUFFER];
static fftw_complex* _fft1d_data[CHANNELS][BUFFER];
static fftw_complex* _fft2d_data[CHANNELS][(ZERO_PAD+STEPS)/2];
static fftw_plan     _fft1d_plan[CHANNELS][BUFFER];
static fftw_plan     _fft2d_plan[CHANNELS][(ZERO_PAD+STEPS)/2];

inline double mag(double real, double imag)
{
  return sqrt((real * real) + (imag * imag));
}

inline double ang(double real, double imag)
{
  return atan(imag / real);
}

void run_init();
void run_daq();
void run_1d_fft();
void run_2d_fft();
void plot_data();

int main()
{

  run_init();

  // control loop
  for (int i = 0; i < 300; i++) {

    run_daq();
    #if FFT==1
    run_1d_fft();
    #endif
    #if FFT==2
    run_1d_fft();
    run_2d_fft();
    #endif
    plot_data();
    //std::thread t_run_daq(run_daq);
    //std::thread t_run_1d_fft(run_1d_fft);
    //std::thread t_plot_data(plot_data, 1);

    //t_plot_data.join();
    //t_run_1d_fft.join();
    //t_run_daq.join();

  }

  return 0;
}

void run_init()
{
  if (!bcm2835_init()) {
    printf("bcm2835_init failed\n");
  }

  if (!bcm2835_spi_begin()) {
    printf("bcm2835_spi_begin failed\n");
  }

  if (!bcm2835_aux_spi_begin()) {
    printf("bcm2835_spi_begin failed\n");
  }

  // adc spi inits
  bcm2835_spi_setBitOrder(BCM2835_SPI_BIT_ORDER_MSBFIRST);
  bcm2835_spi_setDataMode(BCM2835_SPI_MODE0);
  bcm2835_spi_setClockDivider(BCM2835_SPI_CLOCK_DIVIDER_128);
  bcm2835_spi_chipSelect(BCM2835_SPI_CS1);
  bcm2835_spi_setChipSelectPolarity(BCM2835_SPI_CS1, LOW);
    // Conversion pin set to output
  bcm2835_gpio_fsel(RPI_BPLUS_GPIO_J8_33, BCM2835_GPIO_FSEL_OUTP);

  // dac aux inits
  bcm2835_aux_spi_setClockDivider(BCM2835_SPI_CLOCK_DIVIDER_128);

  // allocating space to fftw vars
  for (int c = 0; c < CHANNELS; c++) 
    for (int b = 0; b < BUFFER; b++) {
      _raw_data[c][b] = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * ZERO_PAD_STEPS);
      _fft1d_data[c][b] = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * ZERO_PAD_STEPS);
    }

  for (int c = 0; c < CHANNELS; c++) 
    for (int s = 0; s < ZERO_PAD_STEPS/2; s++)
      _fft2d_data[c][s] = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * BUFFER);


  for (int c = 0; c < CHANNELS; c++) 
    for (int b = 0; b < BUFFER; b++)
      _fft1d_plan[c][b] = fftw_plan_dft_1d(ZERO_PAD_STEPS,
          (fftw_complex*)_raw_data[c][b], (fftw_complex*)_fft1d_data[c][b],
          FFTW_FORWARD, FFTW_ESTIMATE);

  for (int c = 0; c < CHANNELS; c++) 
    for (int s = 0; s < ZERO_PAD_STEPS/2; s++)
      _fft2d_plan[c][s] = fftw_plan_dft_1d(BUFFER,
          (fftw_complex*)_fft2d_data[c][s], (fftw_complex*)_fft2d_data[c][s],
          FFTW_FORWARD, FFTW_ESTIMATE);

  // hamming window
  for (int n = 0; n < STEPS; n++)
  #if HAMMING_WINDOW==true
    hamming_window[n] = 0.54 - (0.46 * cos(2*M_PI*n / STEPS));
  #else
    hamming_window[n] = 1;
  #endif

  for (int s = 0; s < ZERO_PAD_STEPS/2; s++) fft_calibration[s] = 0;
  for (int s = 0; s < STEPS; s++) radar_calibration[s] = 0;

  #if ZERO_FREE_SPACE==true
  
  #endif
  #if FREE_SPACE_CYCLES>0
  printf("\nCalibrating free space for %d cycles...\n", FREE_SPACE_CYCLES);
  
  int sbins = FREE_SPACE_CYCLES * BUFFER;
  
  run_daq();
  
  for (int c = 0; c < FREE_SPACE_CYCLES; c++) {
  
    printf("%d\n", c+1);
    run_daq();
    run_1d_fft();
  
    for (int b = 0; b < BUFFER; b++)
      for (int s = 0; s < STEPS; s++)
        radar_calibration[s] +=  _raw_data[0][b][s][REAL];

    for (int s = 0; s < STEPS; s++)
      radar_calibration[s] = radar_calibration[s] / sbins;

  }
  std::ofstream f("free_space.dat");

  for (int s = 0; s < STEPS; s++)
    f << radar_calibration[s] << '\n';

  f.close();
  #endif
  #if RADAR_CYCLES>0
  printf("\nCalibrating RADAR for %d cycles...\n", RADAR_CYCLES);

  for (int s = 0; s < STEPS; s++) radar_calibration[s] = 0;
  
  int rbins = RADAR_CYCLES * BUFFER;
  
  run_daq();
  
  for (int c = 0; c < RADAR_CYCLES; c++) {
  
    printf("%d\n", c+1);
    run_daq();
    run_1d_fft();
  
    for (int b = 0; b < BUFFER; b++)
      for (int s = 0; s < STEPS; s++)
        radar_calibration[s] +=  _raw_data[0][b][s][REAL] / rbins;

  }
  #endif
  #if FFT_CYCLES>0
  printf("\nCalibrating fft for %d cycles...\n", FFT_CYCLES);
  
  int fbins = FFT_CYCLES * BUFFER;
  
  run_daq();
  
  for (int c = 0; c < FFT_CYCLES; c++) {
  
    printf("%d\n", c+1);
    run_daq();
    run_1d_fft();
  
    for (int b = 0; b < BUFFER; b++)
      for (int s = 0; s < ZERO_PAD_STEPS/2; s++)
        fft_calibration[s] += mag(_fft1d_data[0][b][s][REAL], _fft1d_data[0][b][s][IMAG]) / fbins;
  }
  #endif
}

void run_daq()
{
  static const struct sched_param priority = {1};
  sched_setscheduler(0, SCHED_FIFO, &priority);

  static const int step_delay((1000000 / STEPS / HZ) - TUNE);
  static const uint16_t step_incr(0xFFFF / STEPS);
  static uint16_t cur_mag(0);

  for (int b=0; b<BUFFER; b++, cur_mag=0) {
    dc_offset[b] = 0;
    for (int s=0; s<STEPS; s++) {
      
      /* DAC OUTPUT WRITE */
        // write current magnitude to dac output
      bcm2835_aux_spi_write(cur_mag);

      /* ADC INPUT READ */
        // prepare for adc trade
          // CNV pin -> HIGH
      bcm2835_gpio_set(RPI_BPLUS_GPIO_J8_33);
          // CNV pin has 40ns high time minimum
          // overhead from calling delay fn loosely covers high time
      bcm2835_delayMicroseconds(1); 
          // CNV pin -> LOW
      bcm2835_gpio_clr(RPI_BPLUS_GPIO_J8_33);

        // adc trade
          // trade softspan data for read data on ADC
      bcm2835_spi_transfernb((char*)_softspan_data, (char*)_read_data, (int)SPI_LEN);
          // convert raw data to voltage value
      #if ZERO_FREE_SPACE==true
      _raw_data[0][b][s][REAL] = (double)((0xFFFF\
          - ((_read_data[0] << 8) | (_read_data[1])))\
          * (5.0f/65536.0f)) - free_space_calibration[s];
      #else
      _raw_data[0][b][s][REAL] = (double)((0xFFFF\
          - ((_read_data[0] << 8) | (_read_data[1])))\
          * (5.0f/65536.0f));
      #endif

      //_raw_data[0][p][s][REAL] = (double)((0xFFFF\
      //    - ((_read_data[3] << 8) | (_read_data[4])))\
      //    * (5.0f/65536.0f));
      
      // Get average of the input voltages to calculate DC offset
      #if DC_OFFSET==true
      dc_offset[b] += _raw_data[0][b][s][REAL]/STEPS;
      #endif
      
        // increment DAC magnitude
      cur_mag += step_incr;
        // delay next DAC write to hit target frequency
      bcm2835_delayMicroseconds(step_delay);
      
      //printf("%02x %02x %02x %02x %02x %02x\n",
      //    _read_data[0], _read_data[1], _read_data[2],
      //    _read_data[3], _read_data[4], _read_data[5]
      //    );

    }
  }
  bcm2835_aux_spi_write(0);
}

  //                               //
 /* Run 1d fft over whole dataset */
//                               //
void run_1d_fft()
{
  for (int c = 0; c < CHANNELS; c++)
    for (int b = 0; b < BUFFER; b++) {
      for (int s = 0; s < STEPS; s++)
        _raw_data[c][b][s][REAL] = (_raw_data[c][b][s][REAL]\
          - dc_offset[b] - radar_calibration[s]) * hamming_window[s];
      for (int s = STEPS; s < ZERO_PAD_STEPS; s++)
        _raw_data[0][b][s][REAL] = 0;
    }

  for (int c = 0; c < CHANNELS; c++)
    for (int b = 0; b < BUFFER; b++) fftw_execute(_fft1d_plan[c][b]);
}

  //                                        //
 /* Run 2d fft over half of 1d fft dataset */
//                                        //
void run_2d_fft()
{
  for (int c = 0; c < CHANNELS; c++)
    for (int s = 0; s < ZERO_PAD_STEPS/2; s++)
      for (int b = 0; b < BUFFER; b++) {
        _fft2d_data[c][s][b][REAL] = _fft1d_data[c][b][s][REAL];
        _fft2d_data[c][s][b][IMAG] = _fft1d_data[c][b][s][IMAG];
      }
  for (int c = 0; c < CHANNELS; c++)
    for (int s = 0; s < ZERO_PAD_STEPS/2; s++) fftw_execute(_fft2d_plan[c][s]);
}
  //                      //
 /* Pipe data to gnuplot */
//                      //
void plot_data()
{

  int intr(0);
  static FILE* gp = popen("gnuplot", "w");
  static bool startup = true;
  if (startup == true) {
      // Change stream buffer behavior from line buffering to
      //  full buffering to increase response time
    char buf[10000];
    setvbuf(gp, buf, _IOFBF, sizeof(buf));
    startup = false;
  }
  #if FFT==0
 /* RAW DATA */
//          //
  fprintf(gp,"plot [0:%d][0:6] '-' w lines\n", ZERO_PAD_STEPS);
  for (int b = 0; b < BUFFER; b++)
    for (int s = 0; s < STEPS; s++, intr++)
      fprintf(gp,"%d %lf \n", intr, _raw_data[0][b][s][REAL]);
  fprintf(gp,"e\n");
  fflush(gp);
  #elif FFT==1 && FFT1_PLOT_TYPE==0
 /* 2D 1D FFT PLOT */
//                //
  fprintf(gp,"set grid\nplot [0:%d][%d:%d] '-' pt 5\n", ZERO_PAD_STEPS/2, FFT1_MAG_MIN, FFT1_MAG_MAX);
  for (int b = 0; b < BUFFER; b++)
    for (int s = 0; s < ZERO_PAD_STEPS/2-1; s++, intr++)
      fprintf(gp,"%d %lf \n", intr,
          mag(_fft1d_data[0][b][s][REAL], _fft1d_data[0][b][s][IMAG]) - fft_calibration[s]);
  fprintf(gp,"e\n");
  fflush(gp);
  #elif FFT==1 && FFT1_PLOT_TYPE==1
 /* 3D 1D FFT PLOT */
//                //
  fprintf(gp,"set view map\nset dgrid3d\nset cbrange[%d:%d]\nsplot [0:%d][0:%d] '-' using 1:2:3 with pm3d interpolate 2,2\n",
      FFT1_MAG_MIN, FFT1_MAG_MAX, ZERO_PAD_STEPS/2,BUFFER-1);
  for (int b = 0; b < BUFFER; b++)
    for (int s = 0; s < ZERO_PAD_STEPS/2; s++)
      fprintf(gp, "%d %d %lf\n", s, b,
          mag(_fft1d_data[0][b][s][REAL], _fft1d_data[0][b][s][IMAG])
          - fft_calibration[s]);
  fprintf(gp,"e\n");
  fflush(gp);
  #elif FFT==2 && FFT2_PLOT_TYPE==0
 /* 2D 2D FFT PLOT */
//                //
  static double mag_sum;
  fprintf(gp,"set grid\nplot [0:%d][%d:%d] '-' using 1:2 with lines\n",
               BUFFER-1, FFT2_SUM_MIN, FFT2_SUM_MAX);
  for (int b = BUFFER/2; b < BUFFER; b++) {
    mag_sum = 0;
    for (int s = 0; s < ZERO_PAD_STEPS/2; s++)
      mag_sum += mag(_fft2d_data[0][s][b][REAL], _fft2d_data[0][s][b][IMAG]);
    fprintf(gp, "%d %lf\n", b - BUFFER/2, mag_sum);
    }
  for (int b = CROSSTALK_SAMPLES; b < BUFFER/2; b++) {
    mag_sum = 0;
    for (int s = 0; s < ZERO_PAD_STEPS/2; s++)
      mag_sum += mag(_fft2d_data[0][s][b][REAL], _fft2d_data[0][s][b][IMAG]);
    fprintf(gp, "%d %lf\n", b + BUFFER/2, mag_sum);
  }
  //fprintf(gp,"set multiplot\nset style data histogram\nplot [0:%d][%d:%d] '-' using 1:2 with lines\n",
  //             CROSSTALK_SAMPLES, FFT2_SUM_MIN, FFT2_CT_MAX);
  //for (int b = 0; b < CROSSTALK_SAMPLES; b++) {
  //  mag_sum = 0;
  //  for (int s = 0; s < ZERO_PAD_STEPS/2; s++)
  //    mag_sum += mag(_fft2d_data[0][s][b][REAL], _fft2d_data[0][s][b][IMAG]);
  //  fprintf(gp, "%d %lf\n", b + BUFFER/2, mag_sum);
  //}
  //fprintf(gp, "unset multiplot\n");
  fprintf(gp,"e\n");
  fflush(gp);
  #elif FFT==2 && FFT2_PLOT_TYPE==1
 /* 3D 2D FFT PLOT */
//                //
  //fprintf(gp,"set view map\nset dgrid3d\nset cbrange[%d:%d]\nsplot [0:%d][0:%d] '-' using 1:2:3 with pm3d\n",
  //fprintf(gp,"set dgrid3d\nset pm3d map interpolate 0,0\nset cbrange[%d:%d]\nsplot [0:%d][0:%d] '-' using 1:2:3\n",
  fprintf(gp,"set dgrid3d\nset pm3d map\nset cbrange[%d:%d]\nsplot [0:%d][0:%d] '-' using 1:2:3\n",
              FFT2_MAG_MIN, FFT2_MAG_MAX, BUFFER-1, ZERO_PAD_STEPS/2);
  for (int s = 0; s < ZERO_PAD_STEPS/2; s++) {
    for (int b = BUFFER/2; b < BUFFER; b++)
      fprintf(gp, "%d %d %lf\n", b - BUFFER/2, s,
          mag(_fft2d_data[0][s][b][REAL], _fft2d_data[0][s][b][IMAG]));
    for (int b = 0; b < BUFFER/2; b++)
      fprintf(gp, "%d %d %lf\n", b + BUFFER/2, s,
          mag(_fft2d_data[0][s][b][REAL], _fft2d_data[0][s][b][IMAG]));
  }
  fprintf(gp,"e\n");
  fflush(gp);
  #endif
}
