#include <H5Epublic.h>
#include <H5Spublic.h>

namespace moab { 

extern "C" {
#if defined(H5E_auto_t_vers) && H5E_auto_t_vers > 1
    typedef herr_t (*HDF5_Error_Func_Type)( hid_t, void* );
#else
    typedef herr_t (*HDF5_Error_Func_Type)( void* );
#endif
}

bool HDF5_can_append_hyperslabs();

} // namespace moab

