submodule (ccpp_io_reader) pio_submodule
use pio_reader, only: pio_reader_t
implicit none
contains
  module procedure create_netcdf_reader_t
    allocate(pio_reader_t :: r)
  end procedure create_netcdf_reader_t
end submodule pio_submodule
