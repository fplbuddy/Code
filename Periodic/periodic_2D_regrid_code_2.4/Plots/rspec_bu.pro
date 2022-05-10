N_start=1
N_end=2000
NS=128
spec_av=fltarr(NS/2+1)
FOR isp=N_start,N_end DO BEGIN
                     ext=     strtrim(string(isp),1)
if (isp le 99 ) then ext= '0'+strtrim(string(isp),1)
if (isp le 9  ) then ext='00'+strtrim(string(isp),1)

kk=findgen(NS/2+1)+1
speck=fltarr(NS/2+1)
specm=fltarr(NS/2+1)
close,1
print,  'spectrum_kk.'+ext+'.txt'
openr,1,'spectrum_kk.'+ext+'.txt'
readf,1,speck
close,1
close,1
print,  'spectrum_bb.'+ext+'.txt'
openr,1,'spectrum_bb.'+ext+'.txt'
readf,1,specm
close,1


spec_av=spec_av+speck
 plot,kk,speck*kk*kk*kk*kk,/ylog,/xlog,yr=[1e-12,1e2]
oplot,kk,specm*kk*kk*kk*kk,thick=2
;plot,kk,speck,/xlog,/ylog,yr=[1e-12,1e2]
;plot,kk,spec,/ylog,yr=[1e-12,1e2]
;oplot,kk,specm,thick=2

wait,0.1
ENDFOR
 plot,kk,spec_av,/ylog,/xlog,yr=[1e-8,1e2],thick=4
oplot,kk,1.0/KK^4


END
