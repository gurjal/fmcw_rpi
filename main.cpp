#include "rpi_daq.h"

int main()
{
  rpi_daq radar(16, 16, 128);

  for (int i = 0; i < 10; i++) {

    radar.run_daq();
    radar.run_fft();
    radar.write_plot_data(1);

  }

  return 0;
}
