#include <delay_ops.h>

int main()
{
  std::cout << " Root dir check " << std::endl;
  Delay_Op test("base_parameters.in","smoothing_store.in","high_wind_cap.in","uk_price.in","uk_trading.in","course_grid.in");
  test.print_params(std::cout);

}