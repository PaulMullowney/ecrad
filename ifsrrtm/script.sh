#!/bin/hash
#declare -a arr=(rrtm_taumol1.intfb.h rrtm_taumol2.intfb.h rrtm_taumol3.intfb.h rrtm_taumol4.intfb.h rrtm_taumol5.intfb.h rrtm_taumol6.intfb.h rrtm_taumol7.intfb.h rrtm_taumol8.intfb.h rrtm_taumol9.intfb.h rrtm_taumol10.intfb.h rrtm_taumol11.intfb.h rrtm_taumol12.intfb.h rrtm_taumol13.intfb.h rrtm_taumol14.intfb.h rrtm_taumol15.intfb.h rrtm_taumol16.intfb.h)
#for file in "${arr[@]}"
#do
#    sed "s/INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: laytrop_min, laytrop_max/INTEGER(KIND=JPIM), OPTIONAL, INTENT(INOUT) :: laytrop_min, laytrop_max/g;" $file > temp
#    mv temp $file
#done

declare -a arr=(srtm_taumol16.F90 srtm_taumol17.F90 srtm_taumol18.F90 srtm_taumol19.F90 srtm_taumol20.F90 srtm_taumol21.F90 srtm_taumol22.F90 srtm_taumol23.F90 srtm_taumol24.F90 srtm_taumol25.F90 srtm_taumol26.F90 srtm_taumol27.F90 srtm_taumol28.F90 srtm_taumol29.F90)

#declare -a arr=(srtm_taumol16.F90)
for file in "${arr[@]}"
do
    sed '0,/ )/s/ )/ , laytrop_min, laytrop_max)/' $file > temp
    sed "s/INTEGER(KIND=JPIM) :: laytrop_min, laytrop_max/INTEGER(KIND=JPIM), OPTIONAL, INTENT(INOUT) :: laytrop_min, laytrop_max/g;" temp > temp2
    mv temp2 $file
    rm temp
done
