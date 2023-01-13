;************************************************
begin

ni = 4500
nj = 3297
ncat = 5
nilyr = 4
nslyr = 1

; variables to be written to new restart file

  layers = ispan(1,nilyr,1)
  clayers = sprinti("%0.3i",layers)
  ;print(clayers)
    vars2d = (/      "uvel",       "vvel", "scale_factor",     "coszen" \
              ,     "swvdr",      "swvdf",        "swidr",      "swidf" \
              ,  "strocnxT",   "strocnyT" \
              , "stressp_1",  "stressp_2",    "stressp_3",  "stressp_4" \
              , "stressm_1",  "stressm_2",    "stressm_3",  "stressm_4" \
              ,"stress12_1", "stress12_2",   "stress12_3", "stress12_4" \
              ,"iceumask"/)

  ncatvars = (/"aicen", "vicen", "vsnon", "Tsfcn"/)
  nlyrvars = (/ "qice",  "sice"/)

; set up the output netcdf file

  setfileoption("nc","format","LargeFile")
  foutname="test.nc"
  system("/bin/rm -f " + foutname)    ; remove if exists
  fout = addfile(foutname,"c")

  ; explicitly declare file definition mode. Improve efficiency.
  setfileoption(fout,"DefineMode",True)

  ; create global attributes of the file
  fAtt               = True            ; assign file attributes
  fAtt@creation_date = systemfunc ("date")
  fileattdef( fout, fAtt )           ; copy file attributes

  ; predefine the coordinate variables and their dimensionality
  dimNames = (/"ncat",  "nj",  "ni"/)
  dimSizes = (/ ncat ,    nj,    ni/)
  dimUnlim = (/False , False, False/)
  filedimdef(fout,dimNames,dimSizes,dimUnlim)

  ; predefine the the dimensionality of the variables to be written out
  ;filevardef(fout, "ncat", typeof(ncat), "ncat")
  ;filevardef(fout,   "nj",  typeof(nj),    "nj")
  ;filevardef(fout,   "ni",  typeof(ni),    "ni")

  do nv = 0,dimsizes(ncatvars)-1
    varname = ncatvars(nv)
    odims = (/"ncat", "nj", "ni"/)
    print("creating variable "+varname+" in file")
    filevardef(fout, varname, "double", odims)
    delete(odims)
  end do

  ; predefine variables
  do nv = 0,dimsizes(vars2d)-1
    varname = vars2d(nv)
    odims = (/"nj", "ni"/)
    print("creating variable "+varname+" in file")
    filevardef(fout, varname, "double", odims)
    delete(odims)
  end do

  do nv = 0,dimsizes(nlyrvars)-1
    do nl = 0,nilyr-1
      varname = nlyrvars(nv)+clayers(nl)
      odims = (/"ncat", "nj", "ni"/)
      print("creating variable "+varname+" in file")
      filevardef(fout, varname, "double", odims)
      delete(odims)
    end do
  end do

  ; last two
  varname = "qsno001"
  odims = (/"ncat", "nj", "ni"/)
  print("creating variable "+varname+" in file")
  filevardef(fout, varname, "double", odims)
  delete(odims)

  varname = "frz_onset"
  odims = (/"nj", "ni"/)
  print("creating variable "+varname+" in file")
  filevardef(fout, varname, "double", odims)
  delete(odims)

;if(1 .eq. 0)then

; define binary file and variables to be read/written
fbinname="rtofs_glo.t00z.n-06.restart_cice_converted"
; category variables
aicen = new((/ncat,nj,ni/),double)
vicen = new((/ncat,nj,ni/),double)
vsnon = new((/ncat,nj,ni/),double)
Tsfcn = new((/ncat,nj,ni/),double)

var2d = new((/nj,ni/),double)

setfileoption("bin","WriteByteOrder","BigEndian")

rec=0
do k = 0,ncat-1
  var2d = fbindirread(fbinname,rec,(/nj,ni/),"double")
  print("rec = "+rec+" min val "+min(var2d)+" max val "+max(var2d))
  aicen(k,:,:) = var2d(:,:)
  rec = rec+1

  var2d = fbindirread(fbinname,rec,(/nj,ni/),"double")
  print("rec = "+rec+" min val "+min(var2d)+" max val "+max(var2d))
  vicen(k,:,:) = var2d(:,:)
  rec = rec+1

  var2d = fbindirread(fbinname,rec,(/nj,ni/),"double")
  print("rec = "+rec+" min val "+min(var2d)+" max val "+max(var2d))
  vsnon(k,:,:) = var2d(:,:)
  rec = rec+1

  var2d = fbindirread(fbinname,rec,(/nj,ni/),"double")
  print("rec = "+rec+" min val "+min(var2d)+" max val "+max(var2d))
  Tsfcn(k,:,:) = var2d(:,:)
  rec = rec+1
end do

print ("writing aicen")
fout->aicen = (/aicen/)
;end if

end
