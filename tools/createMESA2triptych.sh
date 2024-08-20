#!/bin/sh

# This script creates a directory to launch a MESA simulation with 
# extra output required for triptych to start a stellar collision.

# First create the directory by copying the defaults from MESA
cp -r $MESA_DIR/star/work ./MESA2triptych

# Enter it
cd ./MESA2triptych

# Add the file profile_columns.list to add additional output to the
# MESA simulation.
cat > profile_columns.list <<-__EOF__
! profile_columns.list -- determines the contents of star model profiles
! you can use a non-standard version by setting profile_columns_file in your inlist

! units are cgs unless otherwise noted.

   zone       ! numbers start with 1 at the surface
   mass       ! m/Msun. mass coordinate of outer boundary of cell.
   logR       ! log10(radius/Rsun) at outer boundary of zone
   logT       ! log10(temperature) at center of zone
   logRho     ! log10(density) at center of zone
   logP       ! log10(pressure) at center of zone
   x_mass_fraction_H
   y_mass_fraction_He
   z_mass_fraction_metals
   h1
   he3
   he4
   c12
   c13
   n14
   n15
   o16
   o17
   o18
   ne20
   mg24
   q

! everything below this line is deactivated
__EOF__

# Finally, add a line to the "inlist" file so that it uses profile_columns.list
sed -i '/&star_job/,/\/ ! end of star_job namelist/ s/\(extra_star_job_inlist_name(1) = '\''inlist_project'\''\)/\1\n    profile_columns_file = '\''profile_columns.list'\''/' inlist

