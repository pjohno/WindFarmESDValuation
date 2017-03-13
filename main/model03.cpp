#include "delay_ops.h"
#include <iostream>

main()
{

  Delay_Op test("base_parameters.in","smoothing_store.in","high_wind_cap.in","uk_price.in","uk_trading.in","course_grid.in");

  test.perpetual_pricer(true,true);

  test.print_params(std::cout);

}