#include <fstream>
#include <vector>
#include <thread>
#include <cmath>
#include <bcm2835.h>
#include <fftw3.h>

#define REAL 0
#define IMAG 1

#define SPI_LEN 6

inline double mag(double real, double imag)
{
  return sqrt((real * real) + (imag * imag));
}

inline double ang(double real, double imag)
{
  return atan(imag / real);
}

class rpi_daq
{
public:
  rpi_daq(int hz, int periods, int steps, int tune = 0);
  ~rpi_daq();
  void run_daq();
  void run_fft();
  void write_plot_data(int n);

private:
  int _hz;
  int _steps;
  int _periods;
  int _tune;
    
  // ADC data for 2 channels
  char _softSpan_data[SPI_LEN] = { 0x00, 0x90, 0x00, 0x00, 0x00, 0x00 };
  char _send_data[SPI_LEN]     = { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };
  char _read_data[SPI_LEN]     = { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };
  
  // FFTW data structures
  fftw_complex ** _ch_A_raw;
  fftw_complex ** _ch_B_raw;
  fftw_complex ** _ch_A_fft;
  fftw_complex ** _ch_B_fft;
  fftw_plan    * _ch_A_plan;
  fftw_plan    * _ch_B_plan;
};

rpi_daq::rpi_daq(int hz, int periods, int steps, int tune)
  : _hz(hz), _steps(steps), _periods(periods), _tune(tune),
    _ch_A_raw(new fftw_complex * [periods]),
    _ch_B_raw(new fftw_complex * [periods]),
    _ch_A_fft(new fftw_complex * [periods]),
    _ch_B_fft(new fftw_complex * [periods]),
    _ch_A_plan(new fftw_plan[periods]),
    _ch_B_plan(new fftw_plan[periods])
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
  for (int p = 0; p < _periods; p++) {
    _ch_A_raw[p] = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * _steps);
    _ch_B_raw[p] = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * _steps);
    _ch_A_fft[p] = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * _steps);
    _ch_B_fft[p] = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * _steps);
  }

  for (int p = 0; p < _periods; p++) {
    _ch_A_plan[p] = fftw_plan_dft_1d(_steps, (fftw_complex*)_ch_A_raw[p],
        (fftw_complex*)_ch_A_fft[p], FFTW_FORWARD, FFTW_ESTIMATE);
    _ch_B_plan[p] = fftw_plan_dft_1d(_steps, (fftw_complex*)_ch_B_raw[p],
        (fftw_complex*)_ch_B_fft[p], FFTW_FORWARD, FFTW_ESTIMATE);
  }
}
    
rpi_daq::~rpi_daq()
{
  delete [] _ch_A_raw;
  delete [] _ch_B_raw;
  delete [] _ch_A_fft;
  delete [] _ch_B_fft;
  delete [] _ch_A_plan;
  delete [] _ch_B_plan;
}

// Write and read DAQ for number of PERIODS
//
void rpi_daq::run_daq() {

  static const int step_delay((1000000 / _steps / _hz) - _tune);
  static const uint16_t step_incr(0xFFFF / _steps);
  static uint16_t cur_mag(0);

  for (int p=0; p<_periods; p++, cur_mag=0) {
    for (int s=0; s<_steps; s++) {
      
      /* DAC OUTPUT WRITE */
        // write current magnitude to dac output
      bcm2835_aux_spi_write(cur_mag);

      /* ADC INPUT READ */
        // prepare for adc trade
          // CNV pin -> HIGH
      bcm2835_gpio_set(RPI_BPLUS_GPIO_J8_33);
          // CNV pin has 40ns high time minimum
          // overhead from calling delay fn loosely covers high time
      bcm2835_delayMicroseconds(0); 
          // CNV pin -> LOW
      bcm2835_gpio_clr(RPI_BPLUS_GPIO_J8_33);

        // adc trade
          // trade softspan data for read data on ADC
      bcm2835_spi_transfernb((char*)_softSpan_data, (char*)_read_data, (int)SPI_LEN);
          // convert raw data to voltage value
      _ch_A_raw[p][s][REAL] = (double)((0xFFFF - ((_read_data[0] << 8) | (_read_data[1]))) * (5.0f/65536.0f));
      _ch_B_raw[p][s][REAL] = (double)((0xFFFF - ((_read_data[3] << 8) | (_read_data[4]))) * (5.0f/65536.0f));
      
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

// Run fft on each period individually
//
void rpi_daq::run_fft() {
  for (int p = 0; p < _periods; p++) {
    fftw_execute(_ch_A_plan[p]);
    fftw_execute(_ch_B_plan[p]);
  }
}

// Write data to file to be plotted by Gnuplot
// Args:
//  0 -> raw adc data
//  1 -> fft on raw adc data
//
void rpi_daq::write_plot_data(int n) {

  int intr(1);
  std::ofstream f("data.dat");

  switch(n)
  {
    case 0:
      for (int p = 0; p < _periods; p++)
        for (int s = 0; s < _steps; s++, intr++)
          f << intr << "\t"
            << _ch_A_raw[p][s][REAL] << "\t"
            << _ch_B_raw[p][s][REAL] << "\n"; 
      f.close();
      break;

    case 1:
      for (int s = 0; s < _steps; s++, intr++)
        f << intr << "\t"
          << _ch_A_fft[1][s][REAL] << "\t"
          << _ch_B_fft[1][s][REAL] << "\n"; 
      f.close();
      break;
  }
}
