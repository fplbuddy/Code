N_start=1
N_end=2000
NS=128
spec_av=fltarr(NS/2+1)
FOR isp=N_start,N_end DO BEGIN
                     ext=     strtrim(string(isp),1)
if (isp le 99 ) then ext= '0'+strtrim(string(isp),1)
if (isp le 9  ) then ext='00'+strtrim(string(isp),1)

kk=findgen(NS/2+1)+1
spec1=fltarr(NS/2+1)
spec2=fltarr(NS/2+1)
spec3=fltarr(NS/2+1)
spec4=fltarr(NS/2+1)
spec5=fltarr(NS/2+1)

close,1
print,  'spectrum_m1.'+ext+'.txt'
openr,1,'spectrum_m1.'+ext+'.txt'
readf,1,spec1
close,1
print,  'spectrum_m2.'+ext+'.txt'
openr,1,'spectrum_m2.'+ext+'.txt'
readf,1,spec2
close,1
print,  'spectrum_m3.'+ext+'.txt'
openr,1,'spectrum_m3.'+ext+'.txt'
readf,1,spec3
close,1
print,  'spectrum_m4.'+ext+'.txt'
openr,1,'spectrum_m4.'+ext+'.txt'
readf,1,spec4
close,1
print,  'spectrum_m5.'+ext+'.txt'
openr,1,'spectrum_m5.'+ext+'.txt'
readf,1,spec5
close,1

 plot,kk,spec1,/ylog,/xlog,yr=[1e-6,1e6]
oplot,kk,spec2,thick=2
oplot,kk,spec3,thick=2
oplot,kk,spec4,thick=2
oplot,kk,spec5,thick=2

wait,0.1
ENDFOR
 plot,kk,spec_av,/ylog,/xlog,yr=[1e-8,1e2],thick=4
oplot,kk,1.0/KK^4


END
