#include <fstream>
#include <iostream>
#include <vector>
#include <thread>
#include <cmath>
#include <bcm2835.h>
#include <fftw3.h>
#include <unistd.h>
#include <sched.h>

#define REAL 0
#define IMAG 1

#define SPI_LEN 6

#define HZ 16
#define STEPS 128
#define PERIODS 16
#define TUNE 0

char _softSpan_data[SPI_LEN] = { 0x00, 0x90, 0x00, 0x00, 0x00, 0x00 };
char _send_data[SPI_LEN]     = { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };
char _read_data[SPI_LEN]     = { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };
fftw_complex**  _ch_A_raw_1(new fftw_complex* [PERIODS]);
fftw_complex**  _ch_B_raw_1(new fftw_complex* [PERIODS]);
fftw_complex**  _ch_A_fft_1(new fftw_complex* [PERIODS]);
fftw_complex**  _ch_B_fft_1(new fftw_complex* [PERIODS]);
fftw_plan*      _ch_A_plan_1(new fftw_plan[PERIODS]);
fftw_plan*      _ch_B_plan_1(new fftw_plan[PERIODS]);
fftw_complex**  _ch_A_raw_2(new fftw_complex* [STEPS]);
fftw_complex**  _ch_B_raw_2(new fftw_complex* [STEPS]);
fftw_complex**  _ch_A_fft_2(new fftw_complex* [STEPS]);
fftw_complex**  _ch_B_fft_2(new fftw_complex* [STEPS]);
fftw_plan*      _ch_A_plan_2(new fftw_plan[STEPS]);
fftw_plan*      _ch_B_plan_2(new fftw_plan[STEPS]);

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
void write_plot_data(int n, int _p = 1);

int main()
{
  static const struct sched_param priority = {1};
  sched_setscheduler(0, SCHED_FIFO, &priority);

  run_init();
  // control
  for (int i = 0; i < 16*100; i++) {
  //while(1) {
  //  if (std::cin.get() == 'q')
  //    break;

    run_daq();
    write_plot_data(0,16);
    //std::thread t_daq(run_daq);
    //std::thread t_plot(write_plot_data, 0, 1);
    //t_plot.join();
    //t_daq.join();

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
  for (int p = 0; p < PERIODS; p++) {
    _ch_A_raw_1[p] = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * STEPS);
    _ch_A_fft_1[p] = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * STEPS);
    _ch_B_raw_1[p] = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * STEPS);
    _ch_B_fft_1[p] = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * STEPS);
  }

  for (int p = 0; p < PERIODS; p++) {
    _ch_A_plan_1[p] = fftw_plan_dft_1d(STEPS,
        (fftw_complex*)_ch_A_raw_1[p], (fftw_complex*)_ch_A_fft_1[p],
        FFTW_FORWARD, FFTW_ESTIMATE);
    _ch_B_plan_1[p] = fftw_plan_dft_1d(STEPS,
        (fftw_complex*)_ch_B_raw_1[p], (fftw_complex*)_ch_B_fft_1[p],
        FFTW_FORWARD, FFTW_ESTIMATE);
  }

  for (int s = 0; s < STEPS; s++) {
    _ch_A_raw_2[s] = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * PERIODS);
    _ch_A_fft_2[s] = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * PERIODS);
    _ch_B_raw_2[s] = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * PERIODS);
    _ch_B_fft_2[s] = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * PERIODS);
  }

  for (int s = 0; s < STEPS; s++) {
    _ch_A_plan_2[s] = fftw_plan_dft_1d(PERIODS,
        (fftw_complex*)_ch_A_raw_2[s], (fftw_complex*)_ch_A_fft_2[s],
        FFTW_FORWARD, FFTW_ESTIMATE);
    _ch_B_plan_2[s] = fftw_plan_dft_1d(PERIODS,
        (fftw_complex*)_ch_B_raw_2[s], (fftw_complex*)_ch_B_fft_2[s],
        FFTW_FORWARD, FFTW_ESTIMATE);
  }
}

void run_daq()
{

  static const int step_delay((1000000 / STEPS / HZ) - TUNE);
  static const uint16_t step_incr(0xFFFF / STEPS);
  static uint16_t cur_mag(0);

  for (int p=0; p<PERIODS; p++, cur_mag=0) {
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
      bcm2835_spi_transfernb((char*)_softSpan_data, (char*)_read_data, (int)SPI_LEN);
          // convert raw data to voltage value
      _ch_A_raw_1[p][s][REAL] = (double)((0xFFFF\
          - ((_read_data[0] << 8) | (_read_data[1])))\
          * (5.0f/65536.0f));
      _ch_B_raw_1[p][s][REAL] = (double)((0xFFFF\
          - ((_read_data[3] << 8) | (_read_data[4])))\
          * (5.0f/65536.0f));
      
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

// Run 1d fft over whole dataset
//
void run_1d_fft()
{
  for (int p = 0; p < PERIODS; p++) {
    fftw_execute(_ch_A_plan_1[p]);
    fftw_execute(_ch_B_plan_1[p]);
  }
}

// Run 2d fft over whole dataset
//
void run_2d_fft()
{
  for (int p = 0; p < PERIODS; p++) {
    fftw_execute(_ch_A_plan_1[p]);
    fftw_execute(_ch_B_plan_1[p]);
  }

  for (int p = 0; p < PERIODS; p++) {
    for (int s = 0; s < STEPS; s++) {
      _ch_A_raw_2[s][p][REAL] = _ch_A_fft_1[p][s][REAL];
      _ch_A_raw_2[s][p][IMAG] = _ch_A_fft_1[p][s][IMAG];
    }
  }

  for (int s = 0; s < STEPS; s++) {
    fftw_execute(_ch_A_plan_2[s]);
    fftw_execute(_ch_B_plan_2[s]);
  }
}

// Write data to file to be plotted by Gnuplot
// Args:
//  0 -> raw adc data
//  1 -> 1d fft data
//  2 -> 2d fft data
//
void write_plot_data(int n, int _p)
{
  int intr(0);
  static FILE* gp = popen("gnuplot", "w");

  switch(n)
  {
    case 0:
      fprintf(gp,"plot [0:%d][0:6] '-' w lines\n", _p*STEPS);
      for (int p = 0; p < _p; p++)
        for (int s = 0; s < STEPS; s++, intr++)
          fprintf(gp,"%d %lf \n", intr, _ch_A_raw_1[p][s][REAL]);
      fprintf(gp,"e\n");
      fflush(gp);
      break;

    case 1:
      fprintf(gp,"plot '-' w lines\n");
      for (int p = 0; p < _p; p++)
        for (int s = 0; s < STEPS; s++, intr++)
          fprintf(gp,"%d %lf \n", intr,
              mag(_ch_A_fft_1[p][s][REAL], _ch_A_fft_1[p][s][IMAG]));
      fprintf(gp,"e\n");
      fflush(gp);
      break;

    //case 2:
      //for (int s = 0; s < _p; s++)
      //  for (int p = 0; p < _periods; p++, intr++)
      //    f << intr << "\t"
      //      << mag(_ch_A_fft_2[s][p][REAL], _ch_A_fft_2[s][p][IMAG]) << "\n";
      //      //<< mag(_ch_B_fft_2[s][p][REAL], _ch_A_fft_2[s][p][IMAG]) << "\n"; 
      //f.close();
      //break;
  }
}
