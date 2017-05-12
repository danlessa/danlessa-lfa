;colorset
;device,decomposed=0
;modct,33,min=0,max=1,/paper

;;--read in SWFLUXANAL from SIRS system--;;
;;--for RACORO time period to use with VARANAL--;;

swflux_files=file_search('arm/senae1/sgpradflux/sgpradflux1longC1.*')
nf=n_elements(swflux_files)

sirs=fltarr(20,1)
;sirs=[0]month,[1]day,[2]hour,[3]min,[4]csza,[5]cldfrac,[6]globaldn,[7]diffusedn,[8]directdn,
;[9]diffglobalratio,[10]dir+diff,[11]globaldn_clearskyfit,[12]directdn_clearskyfit
;[13]surfalbmeasured,[14]surfalbedoclrskyfit,[15]codsirs,[16]b,[17]b1,[18]b2,[19]year
for i=0,nf-1 do begin
    cdfid=ncdf_open(swflux_files(i))
    vn=ncdf_vardir(cdfid)

    varid=ncdf_varid(cdfid,vn(0)) &  ncdf_varget,cdfid,varid,basetime
    varid=ncdf_varid(cdfid,vn(1)) &  ncdf_varget,cdfid,varid,timeoffset

    bt=double(julday(1,1,1970,0)+(basetime/86400.))
    st=(timeoffset(*)/86400.)+bt
    caldat,st,mm,dd,yy,hh,mn,ss

    varid=ncdf_varid(cdfid,vn(4)) &  ncdf_varget,cdfid,varid,glbdn ;global (total) measured
    varid=ncdf_varid(cdfid,vn(6)) &  ncdf_varget,cdfid,varid,glbdnclr ;global (total) clear-sky fit
    varid=ncdf_varid(cdfid,vn(10)) &  ncdf_varget,cdfid,varid,swup ;upwelling sw - new variable
    varid=ncdf_varid(cdfid,vn(12)) &  ncdf_varget,cdfid,varid,swupclr ;upwelling sw clear-sky fit - new variable
    varid=ncdf_varid(cdfid,vn(16)) &  ncdf_varget,cdfid,varid,difdn ;diffuse measured    
    varid=ncdf_varid(cdfid,vn(19)) &  ncdf_varget,cdfid,varid,dirdn ;direct measured
    varid=ncdf_varid(cdfid,vn(21)) &  ncdf_varget,cdfid,varid,dirdnclr ;direct clear-sky fit 
    varid=ncdf_varid(cdfid,vn(24)) &  ncdf_varget,cdfid,varid,cldfrac
    varid=ncdf_varid(cdfid,vn(25)) &  ncdf_varget,cdfid,varid,codsirs ; cod - new variable
    varid=ncdf_varid(cdfid,vn(29)) &  ncdf_varget,cdfid,varid,csza
    varid=ncdf_varid(cdfid,vn(30)) &  ncdf_varget,cdfid,varid,cldtransm ; cloud transmissivity SW - new variable    
    
    varid=ncdf_varid(cdfid,vn(5)) &  ncdf_varget,cdfid,varid,qcglbdnsw ;qc global dn sw
    varid=ncdf_varid(cdfid,vn(11)) &  ncdf_varget,cdfid,varid,qcglbupsw ;qc global up sw
    varid=ncdf_varid(cdfid,vn(17)) &  ncdf_varget,cdfid,varid,qcdiffdnsw ;qc diffuse dn sw
    varid=ncdf_varid(cdfid,vn(20)) &  ncdf_varget,cdfid,varid,qcdirdnsw ;qc diffuse up sw
    
    difratiodn=difdn/glbdn
    sumdn=dirdn+difdn
    salbobs=swup/glbdn
    salbfit=swupclr/glbdnclr
    
    ;;--calculate Xie and Liu 2013 (ERL) cloud albedo and cloud fraction parameters--;;

    b1=(glbdnclr-glbdn)/(glbdnclr-(swup*(cldtransm*cldtransm)))
    b2= (dirdnclr-dirdn)/dirdnclr
    b=b1/b2
    
    ind=where(difratiodn ge 0 and difratiodn le 1 and salbobs ge 0 and salbobs le 1 and b2 gt 0 and b2 le 1 and b gt 0 and b le 1 and qcglbdnsw eq 0 and qcglbupsw eq 0 and qcdiffdnsw eq 0 and qcdirdnsw eq 0, count)
    
    if (count gt 0) then begin
        buf=fltarr(20,n_elements(hh[ind]))
        buf(0,*)=mm[ind]
        buf(1,*)=dd[ind]
        buf(2,*)=hh[ind]
        buf(3,*)=mn[ind]
        buf(4,*)=csza[ind]
        buf(5,*)=cldfrac[ind]
        buf(6,*)=glbdn[ind]
        buf(7,*)=difdn[ind]
        buf(8,*)=dirdn[ind]
        buf(9,*)=difratiodn[ind]
        buf(10,*)=sumdn[ind]
        buf(11,*)=glbdnclr[ind]
        buf(12,*)=dirdnclr[ind]
        buf(13,*)=salbobs[ind]
        buf(14,*)=salbfit[ind]
        buf(15,*)=codsirs[ind]
        buf(16,*)=b[ind]
        buf(17,*)=b1[ind]
        buf(18,*)=b2[ind]
        buf(19,*)=yy[ind]
        sirs=[[sirs],[buf]]
    endif
    ncdf_close,cdfid
    
endfor
sirs=sirs(*,1:*)
dex=where(sirs(5,*) lt 0,cnt) & if (cnt ge 1) then sirs(5,dex)=!Values.F_NaN
dex=where(sirs(6,*) le 0,cnt) & if (cnt ge 1) then sirs(6,dex)=!Values.F_NaN
dex=where(sirs(7,*) le 0,cnt) & if (cnt ge 1) then sirs(7,dex)=!Values.F_NaN
dex=where(sirs(8,*) le 0,cnt) & if (cnt ge 1) then sirs(8,dex)=!Values.F_NaN
dex=where(sirs(9,*) lt 0,cnt) & if (cnt ge 1) then sirs(9,dex)=!Values.F_NaN
dex=where(sirs(10,*) le 0,cnt) & if (cnt ge 1) then sirs(10,dex)=!Values.F_NaN
dex=where(sirs(11,*) le 0,cnt) & if (cnt ge 1) then sirs(11,dex)=!Values.F_NaN
dex=where(sirs(12,*) le 0,cnt) & if (cnt ge 1) then sirs(12,dex)=!Values.F_NaN
dex=where(sirs(13,*) le 0,cnt) & if (cnt ge 1) then sirs(13,dex)=!Values.F_NaN
dex=where(sirs(14,*) le 0,cnt) & if (cnt ge 1) then sirs(14,dex)=!Values.F_NaN
dex=where(sirs(15,*) lt 0,cnt) & if (cnt ge 1) then sirs(15,dex)=!Values.F_NaN

save,sirs,file='swflux_sirs1min.sav'

g=0.87 ;assumed asymmetry parameter (Xie&Liu 2013 from Hu&Stamnes 1993)   
xb=sirs(16,*) & xb1=sirs(17,*) & xb2=sirs(18,*) & cossza=sirs(4,*)
a_r=fltarr(n_elements(xb)) & a_r(*)=!Values.F_NaN
ma_r=fltarr(n_elements(xb)) & ma_r(*)=!Values.F_NaN
f=fltarr(n_elements(xb)) & f(*)=!Values.F_NaN
t=fltarr(n_elements(xb)) & t(*)=!Values.F_NaN
for j=0,n_elements(xb)-1 do begin
    nb1=xb1(j) & nb2=xb2(j) & nb=xb(j)
    if (finite(nb) ne 0 and finite(nb1) ne 0 and finite(nb2) ne 0) then begin
      if (nb lt 0.07872) then a_r(j)=0
      if (nb ge 0.07872 and nb lt 0.11442) then a_r(j)=1-31.1648*nb+sqrt(((31.1648*nb)^2)-(49.6255*nb))
      if (nb ge 0.11442 and nb le 0.185) then a_r(j)=(((2.61224*nb1)-nb2)+sqrt((24.2004*(nb1^2))-(9.0098*nb1*nb2)+(nb2^2)))/((18.3622*nb1)-(4*nb2))
      if (nb gt 0.185 and nb le 0.23792) then a_r(j)=(0.89412*nb)+0.02519
      if (nb gt 0.23792 and nb le 1.0) then a_r(j)=nb
    endif
    ;;modify cld albedo for neglect of absorption
    ma_r(j)=a_r(j)/(1.0537+(0.0788*cossza(j)))
    if a_r(j) ne 0 then f(j)=xb1(j)/ma_r(j) ;cloud fraction
    t(j)=(2*ma_r(j)*cossza(j))/((1-ma_r(j))*(1-g))  
endfor

acrf=sirs(6,*)-sirs(11,*) ;abs cloud radiative forcing (total dn - clear dn)
rcrf=1-(sirs(6,*)/sirs(11,*)) ;relative cloud radiative forcing (Liu et al.2011)


;crf=[0]month,[1]day,[2]hour,[3]min,[4]csza,[5]cldfrac-sirs,[6]globaldn,[7]diffusedn,[8]directdn,
;[9]crf,[10]rcrf,[11]cloud albedo,[12]cldfrac-Xie&Liu,[13]cloudoptical depth-Xie&Liu,[14]year

crf=fltarr(15,n_elements(sirs(0,*)))
crf(0:8,*)=sirs(0:8,*)
crf(9,*)=acrf
crf(10,*)=rcrf
crf(11,*)=ma_r
crf(12,*)=f
crf(13,*)=t
crf(14,*)=sirs(19,*)

end
