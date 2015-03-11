
## 20150311
* When using snowbands, band area weighted elevation cell elevation must equal 
  that of elevation parameter file to within +/-0.1 m 
* driver/SFLXALL_SRC.f>ALBCPU(COEF=0.5): adjusts the dynamic snow albedo to be an average of prescribed snow albedo 
  (input as parameter by calendar month) and that of albedo decay function (with snow aging). This was set to 1.0
  and thus fully exlcuded **either the decay or monthly albedos** 
* driver/SFLXALL_SRC.f>REDPRM(SNUPX): by vegetation type, adjusts the snow depth threshold at which maximum snow albedo is allowed.  
  A generally shallower threshold was in place and a deeper set of thresholds was commented out.
  Per Ben's recommendation, the deeper set was re-instituted.  

## 20141209
* Reduced number of output variables for purpose of Integrated Scenarios production.  
  See driver/{OPEN,WRITE}_OUTPUT.f90.ISprod.
* Added default Sept-1:00h snow scraping that can be overrode by control file SCRAPESNOW=0.
* Recent change to driver/READ_INITIAL.f90 so that initial state file is read correctly.
  Stripped out Livneh code that assumed warm restart (neglected frozen soils).
* Recent change to physics/SFLXALL_SRC.f>REDPRM() that isolated/quarantined code block 
  near end that set NROOT=0 which garbled TRANSP() and others.

