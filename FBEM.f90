program main
      
        implicit none
      ! will use name list for key parameters
        real,parameter::Ox=0.21,Ca=410,Rgas=8.314 !Ca in ppm
        real:: Ta=20.0,swc=0.35,I=500.0,RH=0.6,LAI=3.0 !driving force;I: absorbed PAR
        real:: f_Ci,Ci,Tk,alpha_q
        real:: Vm_25,Ea_vm,Gamma_star25,Ea_gamma,Ea_Kc,Ea_Ko,Kc_25,Ko_25
        real:: r_JmVm
        real:: gl,D0,kn
        real:: Reco0,Q10,a1
        real:: Vm,Gamma_star,Kc,Ko,Jm,Jc,Je,A
        real:: es,D,Gs,An,Ac,Reco,NEE
        namelist /main_para/ f_Ci,alpha_q,Vm_25,Ea_vm,Gamma_star25,Ea_gamma,Ea_Kc, &
                             Ea_Ko,Kc_25,Ko_25,r_JmVm,gl,D0,kn,Reco0,Q10,a1

        open(11,file='FBEM_namelist.nml')
        read(11,nml=main_para)
        close(11)
        
        Ci = f_Ci * Ca
        Tk = Ta+273.15
        Vm = Vm_25 * Arrhenius(Ea_vm,Rgas,Tk)
        Gamma_star = Gamma_star25*Arrhenius(Ea_gamma,Rgas,Tk)
        Kc = Kc_25*Arrhenius(Ea_Kc,Rgas,Tk)
        Ko = Ko_25*Arrhenius(Ea_Ko,Rgas,Tk)
        Jm = r_JmVm*Vm

        Jc = Vm*(Ci-Gamma_star)/(Ci+Kc*(1.+Ox/Ko))
        Je = (alpha_q*I*Jm/(sqrt(Jm**2+alpha_q**2*I**2)))*((Ci-Gamma_star)/(4*(Ci+2*Gamma_star)))  
        A = min(Jc,Je)
        
        es = exp(21.382-5347.5/Tk)
        D = 0.1*es*(1-RH)

        Gs = gl*A/((Ca-Gamma_star)*(1+D/D0))
        An = Gs*(Ca-Ci)! top layer canopy photosynthesis
        Ac = An*(1.0-exp(-kn*LAI))/kn

        Reco = Reco0*Q10**(Ta/10)*(swc/(swc+a1)) ! Ta in degree C;a1:moisture coefficient at which respiration is half the maximum
        NEE = Reco-Ac
        
        print*,'Vm=',Vm, 'A=',A
        print*,'Gs=',Gs, 'D=',D
        print*,'Jc=',Jc,'Je=',Je
        print*,'Reco=',Reco,'Ac=',Ac
        print*,'An=',An,'NEE= ',NEE
    contains

!      function Jc(Vmx,Ci,Gamma_star,Kc,Ox,Ko)
!             real:: Vmx,Ci,Gamma_star,Kc,Ox,Ko
!              real:: Jc
!              Jc=Vmx*(Ci-Gamma_star)/(Ci+Kc*(1.+Ox/Ko))
!      end function Jc

      function Arrhenius(Ea,Rgas,Tk)
               real:: Ea,Tk
               real:: Rgas 
               real:: Arrhenius
               Arrhenius=exp((Ea/(Rgas*Tk)*(Tk/298.15-1.)))
      end function Arrhenius


end program main
