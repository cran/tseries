.First.lib <- function (lib, pkg)
{
  library.dynam("tseries", pkg, lib)
  provide(tseries)
  require(ts)
}
