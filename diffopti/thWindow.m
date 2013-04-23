function Ttw = thWindow ( list, params, Rmin, Rtox )
%WINDOW Summary of this function goes here
%   Detailed explanation goes here
    Co = 100;
    Avd = 3.0272;
    Bvd = -2.2772;
    Cvd = 0.0875;
    Dvd = -0.08;
    D01 = 1.27e+4;
    Kr = 0.4;
    dens = 1145;
    hp = 0.2;
    tf = 15000;
    dt = 10;

    li=1;
    if( ismember('Conc',params) ) 
      Co=list(li);
      li=li+1;
    end
    if( ismember('A',params) ) 
      Avd=list(li);
      li=li+1;
    end
    if( ismember('B',params) ) 
      Bvd=list(li);
      li=li+1;
    end
    if( ismember('C',params) ) 
      Cvd=list(li);    
      li=li+1;
    end
    if( ismember('D',params) ) 
      Dvd=list(li);   
      li=li+1;
    end
    if( ismember('D01',params) ) 
      D01=list(li);   
      li=li+1;
    end
    if( ismember('poly partition',params) ) 
      Kr=list(li);   
      li=li+1;
    end
    if( ismember('density',params) ) 
      dens=list(li);   
      li=li+1;
    end
    if( ismember('thickness',params) ) 
      hp=list(li);   
      li=li+1;
    end
    if( ismember('duration',params) ) 
      tf=list(li);   
      li=li+1;
    end
    if( ismember('time step',params) ) 
      dt=list(li);   
      li=li+1;
    end

    M=release(Co, Avd, Bvd, Cvd, Dvd, D01, Kr, dens, hp, tf, dt);
    
    n=M.result.numerical('int2').getReal;
    t1=0:20:1000;
    t2=1000:1000:tf;
    t=[t1 t2];
    i=1;
    while(i<=length(n) && n(i)<Rmin)
        i=i+1;
    end
    tdi=i;
    i=length(n);
    while(i>0 && n(i)<Rmin)
        i=i-1;
    end
    tdd=i;

    if(tdi>length(n))
      tdi=length(n);
    end
    if(tdd<1)
      tdd=1;
    end

    Tinc=t(tdi);
    Tdec=t(tdd);
    Ttw=Tdec-Tinc;
    if(Ttw<=0)
        Ttw=0;    
    else
      i=1;
      while(i<=length(n))
	  if(n(i)>=Rtox)
	      Ttw=0;
	  end
	  i=i+1;
      end

      if(Ttw>0)
	if(tdi > 1)
	    M.result.numerical('int2').set('interp', {strcat('range(',num2str(t(tdi-1)),', 1 ,',num2str(t(tdi)),')')});
	    n1=M.result.numerical('int2').getReal;
	    minI=1;
	    i=1;
	    while(i<=length(n1))
		if(abs(Rmin-n1(minI)) > abs(Rmin-n1(i)))
		    minI=i;
		end
		i=i+1;
	    end
	    Tinc=t(tdi-1)+1*(minI-1);
	end
      
	if(tdd < length(n))
	    M.result.numerical('int2').set('interp', {strcat('range(',num2str(t(tdd)),', 1 ,',num2str(t(tdd+1)),')')});
	    n1=M.result.numerical('int2').getReal;
	    minI=1;
	    i=1;
	    while(i<=length(n1))
		if(abs(Rmin-n1(minI)) > abs(Rmin-n1(i)))
		    minI=i;
		end
		i=i+1;
	    end
	    Tdec=t(tdd)+1*(minI-1);
	end
      
      Ttw=Tdec-Tinc;
    end
  end

end

