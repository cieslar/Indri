#!/bin/bash

#wget -O - "http://www.atnf.csiro.au/research/pulsar/psrcat/proc_form.php?version=1.52&Name=Name&JName=JName&GL=GL&GB=GB&P0=P0&P1=P1&DM=DM&RM=RM&Binary=Binary&Dist=Dist&Dist_DM=Dist_DM&Survey=Survey&OSurvey=OSurvey&R_lum=R_lum&R_lum14=R_lum14&Age=Age&Bsurf=Bsurf&P1_i=P1_i&Age_i=Age_i&Bsurf_i=Bsurf_i&startUserDefined=true&c1_val=&c2_val=&c3_val=&c4_val=&sort_attr=Binary&sort_order=asc&condition=&pulsar_names=&ephemeris=short&coords_unit=raj%2Fdecj&radius=&coords_1=&coords_2=&style=Long+with+errors&no_value=nan&fsize=3&x_axis=&x_scale=linear&y_axis=&y_scale=linear&state=query&table_bottom.x=44&table_bottom.y=27" > ATNF_MOD1.html

wget -O - "http://www.atnf.csiro.au/research/pulsar/psrcat/proc_form.php?version=1.54&Name=Name&JName=JName&GL=GL&GB=GB&P0=P0&P1=P1&DM=DM&RM=RM&Binary=Binary&Dist=Dist&Dist_DM=Dist_DM&Survey=Survey&OSurvey=OSurvey&S400=S400&S1400=S1400&Age=Age&Bsurf=Bsurf&P1_i=P1_i&Age_i=Age_i&Bsurf_i=Bsurf_i&startUserDefined=true&c1_val=&c2_val=&c3_val=&c4_val=&sort_attr=Binary&sort_order=asc&condition=&pulsar_names=&ephemeris=short&coords_unit=raj%2Fdecj&radius=&coords_1=&coords_2=&style=Long+with+errors&no_value=nan&fsize=3&x_axis=&x_scale=linear&y_axis=&y_scale=linear&state=query&table_bottom.x=44&table_bottom.y=27" > ATNF_MOD1.html



sed -n -e "/^[0-9]/p" ATNF_MOD1.html > ATNF_MOD1.txt
sed -i 's/<\/a>&nbsp//g' ATNF_MOD1.txt
sed -i 's/<\/a>//g' ATNF_MOD1.txt
sed -i 's/<[^>]*>//g' ATNF_MOD1.txt
sed -i "s/NAN/nan/g" ATNF_MOD1.txt

awk '{if($20=="nan"){print $0}}' ATNF_MOD1.txt >temp
mv temp ATNF_MOD1.txt
#odrzucam te z ujemna Pdot
awk '{if($11>0.){print $0}}' ATNF_MOD1.txt >temp
mv temp ATNF_MOD1.txt
#Pomin recyklingowane pulsary
awk '{if($34>1e10){print $0}}' ATNF_MOD1.txt >temp
mv temp ATNF_MOD1.txt

if [ -d ../input ]
then
	mv ATNF_MOD1.txt ../input
	rm -fr  ATNF_MOD1.html
fi
