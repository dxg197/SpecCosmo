/*@@ Startup file for thorn SpecGRMHD @@*/

#include "cctk.h"

static const char *rcsid = "$Header: /cactusdevcvs/SpecCosmo/SpecGRMHD/src/Startup.c Exp $";

CCTK_FILEVERSION(SpecCosmo_SpecGRMHD_Startup_c)

int SpecGRMHD_Startup(void);

int SpecGRMHD_Startup(void)
{

   const char *banner = "Early Universe Simulator";

   CCTK_RegisterBanner(banner);

   return 0;
}
