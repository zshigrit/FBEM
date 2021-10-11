program main
      
        implicit none
      ! will use name list for key parameters
        real,parameter::Ox=0.21,Ca=410,Rgas=8.314 !Ca in ppm
        real:: Ta=20.0,swc=0.35,I=1500.0,RH=0.6,LAI=0.0 !driving force;I: absorbed PAR
        real:: Rad
        real:: f_Ci,Ci,Tk,alpha_q
        real:: Vm_25,Ea_vm,Gamma_star25,Ea_gamma,Ea_Kc,Ea_Ko,Kc_25,Ko_25
        real:: r_JmVm
        real:: g1,D0,kn
        real:: Reco0,Q10,a1
        real:: Vm,Gamma_star,Kc,Ko,Jm,Jc,Je,A
        real:: es,D,Gs,An,Ac,Reco,NEE
        character:: header
        integer:: istat3,day,hour
        namelist /main_para/ f_Ci,alpha_q,Vm_25,Ea_vm,Gamma_star25,Ea_gamma,Ea_Kc, &
                             Ea_Ko,Kc_25,Ko_25,r_JmVm,g1,D0,kn,Reco0,Q10,a1

        open(11,file='FBEM_namelist.nml')
        read(11,nml=main_para)
        close(11)
        
        open(12,file='tair_rh_rad_hourly.in')!,status='old',access ='sequential',form='formatted')

        read(12,'(a100)') header

  
        print*,header
        

        close(12)
        close(13)
        close(14)
        !close(15)
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
