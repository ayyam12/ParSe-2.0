c----------------------------------------------------------------
c
c Mann Whitney U Test compares two sets at a time, and only
c for one property.
c
c THE SET WITH THE SMALLEST NUMBER OF ENTRIES HAS TO GO FIRST.
c
c The output is the one-tail p-value.
c
c STW 10/15/2021
c
c----------------------------------------------------------------

      program Find_best_sequence_property_Mann_Whitney_U_test
      implicit none

      integer length,num,i,j,npep,jj,
     &        num_ala,num_cys,num_asp,num_glu,num_phe,
     &        num_gly,num_his,num_ile,num_lys,num_leu,
     &        num_met,num_asn,num_pro,num_gln,num_arg,
     &        num_ser,num_thr,num_val,num_trp,num_tyr,
     &        num_set,num_set2,num_AA_props,best_AA_prop
      character code_inp*11000,code(11000)*3,remark*3
      real pppiia_h,pppiip_h,pppiig_h,pppiic_h,pppiid_h,pppiie_h,
     &     pppiif_h,pppiih_h,pppiii_h,pppiik_h,pppiil_h,pppiim_h,
     &     pppiin_h,pppiiq_h,pppiir_h,pppiis_h,pppiit_h,pppiiv_h,
     &     pppiiw_h,pppiiy_h,net_charge,
     &     sum_ppii,v_exponent,fppii,turn,v_flory,rh,
     &     scaling_exp(1000),AA_prop(1000),fb_hps(1000),
     &     a_prop,c_prop,d_prop,U,Z,rank(1000),
     &     e_prop,f_prop,g_prop,h_prop,
     &     i_prop,k_prop,l_prop,m_prop,
     &     n_prop,p_prop,q_prop,r_prop,
     &     s_prop,t_prop,v_prop,w_prop,
     &     y_prop,tmp_min,tmp_max,
     &     fb_hps_a,fb_hps_p,fb_hps_g,fb_hps_c,fb_hps_d,fb_hps_e,
     &     fb_hps_f,fb_hps_h,fb_hps_i,fb_hps_k,fb_hps_l,fb_hps_m,
     &     fb_hps_n,fb_hps_q,fb_hps_r,fb_hps_s,fb_hps_t,fb_hps_v,
     &     fb_hps_w,fb_hps_y,sum_fb_hps

      double precision z_score,erf,erfc,p_value,min_p_value


c PPII bias measured in peptides by Hilser group
c Prot Sci 2013, vol 22, pgs 405-417, in a table in supplementary information

      pppiia_h=0.37
      pppiic_h=0.25
      pppiid_h=0.30
      pppiie_h=0.42
      pppiif_h=0.17
      pppiig_h=0.13
      pppiih_h=0.20
      pppiii_h=0.39
      pppiik_h=0.56
      pppiil_h=0.24
      pppiim_h=0.36
      pppiin_h=0.27
      pppiip_h=1.00
      pppiiq_h=0.53
      pppiir_h=0.38
      pppiis_h=0.24
      pppiit_h=0.32
      pppiiv_h=0.39
      pppiiw_h=0.25
      pppiiy_h=0.25

c hydrophobicity scale optimized for LLPS using simulation methods (Robert Best, J Phys Chem B (2021))

      fb_hps_a=0.51507
      fb_hps_c=0.46169
      fb_hps_d=0.30525
      fb_hps_e=0.42621
      fb_hps_f=1.17854
      fb_hps_g=1.24153
      fb_hps_h=0.55537
      fb_hps_i=0.83907
      fb_hps_k=0.47106
      fb_hps_l=0.51207
      fb_hps_m=0.64648
      fb_hps_n=0.78447
      fb_hps_p=0.34128
      fb_hps_q=0.29516
      fb_hps_r=0.24025
      fb_hps_s=0.11195
      fb_hps_t=0.27538
      fb_hps_v=0.55645
      fb_hps_w=0.97588
      fb_hps_y=1.04266

c cycle thru the different amino acid properties in the amino acid index.

      num_AA_props=0
      best_AA_prop=0
      min_p_value=100.0
      open (8,file='aaindex1_new.txt',status='old')
80    continue
      read(8,*,err=15,end=800) remark
      if(remark.eq.'I  ') then
      num_AA_props=num_AA_props+1
      read(8,*)a_prop,r_prop,n_prop,d_prop,
     &         c_prop,q_prop,e_prop,g_prop,
     &         h_prop,i_prop
      read(8,*)l_prop,k_prop,m_prop,f_prop,
     &         p_prop,s_prop,t_prop,w_prop,
     &         y_prop,v_prop

c      if (num_AA_props.eq.1) then
c      write(*,'(10f8.3)')a_prop,r_prop,n_prop,d_prop,
c     &         c_prop,q_prop,e_prop,g_prop,
c     &         h_prop,i_prop
c      write(*,'(10f8.3)')l_prop,k_prop,m_prop,f_prop,
c     &         p_prop,s_prop,t_prop,w_prop,
c     &         y_prop,v_prop
c      stop
c      endif
c      goto 80
c      if (num_AA_props.ne.162) goto 80

c Normalize each amino acid scale to have values from 0 to 1
      tmp_min=1000.0
      if (a_prop.lt.tmp_min) tmp_min=a_prop
      if (c_prop.lt.tmp_min) tmp_min=c_prop
      if (d_prop.lt.tmp_min) tmp_min=d_prop
      if (e_prop.lt.tmp_min) tmp_min=e_prop
      if (f_prop.lt.tmp_min) tmp_min=f_prop
      if (g_prop.lt.tmp_min) tmp_min=g_prop
      if (h_prop.lt.tmp_min) tmp_min=h_prop
      if (i_prop.lt.tmp_min) tmp_min=i_prop
      if (k_prop.lt.tmp_min) tmp_min=k_prop
      if (l_prop.lt.tmp_min) tmp_min=l_prop
      if (m_prop.lt.tmp_min) tmp_min=m_prop
      if (n_prop.lt.tmp_min) tmp_min=n_prop
      if (p_prop.lt.tmp_min) tmp_min=p_prop
      if (q_prop.lt.tmp_min) tmp_min=q_prop
      if (r_prop.lt.tmp_min) tmp_min=r_prop
      if (s_prop.lt.tmp_min) tmp_min=s_prop
      if (t_prop.lt.tmp_min) tmp_min=t_prop
      if (v_prop.lt.tmp_min) tmp_min=v_prop
      if (w_prop.lt.tmp_min) tmp_min=w_prop
      if (y_prop.lt.tmp_min) tmp_min=y_prop
      a_prop=a_prop-tmp_min
      c_prop=c_prop-tmp_min
      d_prop=d_prop-tmp_min
      e_prop=e_prop-tmp_min
      f_prop=f_prop-tmp_min
      g_prop=g_prop-tmp_min
      h_prop=h_prop-tmp_min
      i_prop=i_prop-tmp_min
      k_prop=k_prop-tmp_min
      l_prop=l_prop-tmp_min
      m_prop=m_prop-tmp_min
      n_prop=n_prop-tmp_min
      p_prop=p_prop-tmp_min
      q_prop=q_prop-tmp_min
      r_prop=r_prop-tmp_min
      s_prop=s_prop-tmp_min
      t_prop=t_prop-tmp_min
      v_prop=v_prop-tmp_min
      w_prop=w_prop-tmp_min
      y_prop=y_prop-tmp_min
      tmp_max=-10.0
      if (a_prop.gt.tmp_max) tmp_max=a_prop
      if (c_prop.gt.tmp_max) tmp_max=c_prop
      if (d_prop.gt.tmp_max) tmp_max=d_prop
      if (e_prop.gt.tmp_max) tmp_max=e_prop
      if (f_prop.gt.tmp_max) tmp_max=f_prop
      if (g_prop.gt.tmp_max) tmp_max=g_prop
      if (h_prop.gt.tmp_max) tmp_max=h_prop
      if (i_prop.gt.tmp_max) tmp_max=i_prop
      if (k_prop.gt.tmp_max) tmp_max=k_prop
      if (l_prop.gt.tmp_max) tmp_max=l_prop
      if (m_prop.gt.tmp_max) tmp_max=m_prop
      if (n_prop.gt.tmp_max) tmp_max=n_prop
      if (p_prop.gt.tmp_max) tmp_max=p_prop
      if (q_prop.gt.tmp_max) tmp_max=q_prop
      if (r_prop.gt.tmp_max) tmp_max=r_prop
      if (s_prop.gt.tmp_max) tmp_max=s_prop
      if (t_prop.gt.tmp_max) tmp_max=t_prop
      if (v_prop.gt.tmp_max) tmp_max=v_prop
      if (w_prop.gt.tmp_max) tmp_max=w_prop
      if (y_prop.gt.tmp_max) tmp_max=y_prop
      a_prop=a_prop/tmp_max
      c_prop=c_prop/tmp_max
      d_prop=d_prop/tmp_max
      e_prop=e_prop/tmp_max
      f_prop=f_prop/tmp_max
      g_prop=g_prop/tmp_max
      h_prop=h_prop/tmp_max
      i_prop=i_prop/tmp_max
      k_prop=k_prop/tmp_max
      l_prop=l_prop/tmp_max
      m_prop=m_prop/tmp_max
      n_prop=n_prop/tmp_max
      p_prop=p_prop/tmp_max
      q_prop=q_prop/tmp_max
      r_prop=r_prop/tmp_max
      s_prop=s_prop/tmp_max
      t_prop=t_prop/tmp_max
      v_prop=v_prop/tmp_max
      w_prop=w_prop/tmp_max
      y_prop=y_prop/tmp_max

c      write(*,'(10f8.3)')a_prop,r_prop,n_prop,d_prop,
c     &         c_prop,q_prop,e_prop,g_prop,
c     &         h_prop,i_prop
c      write(*,'(10f8.3)')l_prop,k_prop,m_prop,f_prop,
c     &         p_prop,s_prop,t_prop,w_prop,
c     &         y_prop,v_prop
c      stop


c Read input sequences for first set of proteins, using fasta format.

      open (2,file='ID_sequences_133.fasta',status='old',err=70)
      num_set=0
      goto 10
70    write (*,'("Input file for first sequence set not present")')
      stop
22    format (A11000)

10    continue
      read(2,22,err=11,end=200) code_inp
      length=len(code_inp)
      if(code_inp(1:1).eq.'>') goto 10

c  convert any lower case letter to upcase.
    
      do i=1,length
         num=ichar(code_inp(i:i))
         if (num.ge.97.and.num.le.122) code_inp(i:i) = char(num-32)
      enddo

c  Determine sequence length.

      j=0
      do i=1,length
      if (code_inp(i:i).eq.' ') goto 1
         j=j+1
         code(j)=code_inp(i:i)
1     continue
      enddo 
      npep=j

100   read(2,22,err=11,end=200) code_inp
      length=len(code_inp)
      if(code_inp(1:1).eq.'>') goto 200
      if(npep.ge.10000) goto 100

c  Add to sequence length for multiline entries.
    
      do i=1,length
         num=ichar(code_inp(i:i))
         if (num.ge.97.and.num.le.122) code_inp(i:i) = char(num-32)
      enddo

      j=npep
      do i=1,length
      if (code_inp(i:i).eq.' ') goto 2
         j=j+1
         code(j)=code_inp(i:i)
2     continue
      enddo
      npep=j
      
      goto 100

200   continue
      if(npep.gt.10000) goto 10
      num_set=num_set+1

      fb_hps(num_set)=0.0
      scaling_exp(num_set)=0.0
      AA_prop(num_set)=0.0

c calculate v and turn propensity

      num_ala=0
      num_cys=0
      num_asp=0
      num_glu=0
      num_phe=0
      num_gly=0
      num_his=0
      num_ile=0
      num_lys=0
      num_leu=0
      num_met=0
      num_asn=0
      num_pro=0
      num_gln=0
      num_arg=0
      num_ser=0
      num_thr=0
      num_val=0
      num_trp=0
      num_tyr=0

      DO J=1,NPEP
         IF (CODE(J).EQ.'A') THEN
            num_ala=num_ala+1
         endif
         IF (CODE(J).EQ.'C') THEN
            num_cys=num_cys+1
         endif
         IF (CODE(J).EQ.'D') THEN
            num_asp=num_asp+1
         endif
         IF (CODE(J).EQ.'E') THEN
            num_glu=num_glu+1
         endif
         IF (CODE(J).EQ.'F') THEN
            num_phe=num_phe+1
         endif
         IF (CODE(J).EQ.'G') THEN
            num_gly=num_gly+1
         endif
         IF (CODE(J).EQ.'H') THEN
            num_his=num_his+1
         endif
         IF (CODE(J).EQ.'I') THEN
            num_ile=num_ile+1
         endif
         IF (CODE(J).EQ.'K') THEN
            num_lys=num_lys+1
         endif
         IF (CODE(J).EQ.'L') THEN
            num_leu=num_leu+1
         endif
         IF (CODE(J).EQ.'M') THEN
            num_met=num_met+1
         endif
         IF (CODE(J).EQ.'N') THEN
            num_asn=num_asn+1
         endif
         IF (CODE(J).EQ.'P') THEN
            num_pro=num_pro+1
         endif
         IF (CODE(J).EQ.'Q') THEN
            num_gln=num_gln+1
         endif
         IF (CODE(J).EQ.'R') THEN
            num_arg=num_arg+1
         endif
         IF (CODE(J).EQ.'S') THEN
            num_ser=num_ser+1
         endif
         IF (CODE(J).EQ.'T') THEN
            num_thr=num_thr+1
         endif
         IF (CODE(J).EQ.'V') THEN
            num_val=num_val+1
         endif
         IF (CODE(J).EQ.'W') THEN
            num_trp=num_trp+1
         endif
         IF (CODE(J).EQ.'Y') THEN
            num_tyr=num_tyr+1
         endif
      enddo

      net_charge=abs(num_asp+num_glu
     &                  -num_lys-num_arg)

      sum_ppii=0.0

      sum_ppii=num_ala*pppiia_h+num_cys*pppiic_h
     &        +num_asp*pppiid_h+num_glu*pppiie_h
     &        +num_phe*pppiif_h+num_gly*pppiig_h
     &        +num_his*pppiih_h+num_ile*pppiii_h
     &        +num_lys*pppiik_h+num_leu*pppiil_h
     &        +num_met*pppiim_h+num_asn*pppiin_h
     &        +num_pro*pppiip_h+num_gln*pppiiq_h
     &        +num_arg*pppiir_h+num_ser*pppiis_h
     &        +num_thr*pppiit_h+num_val*pppiiv_h
     &        +num_trp*pppiiw_h+num_tyr*pppiiy_h

      fppii=sum_ppii/real(npep)
      v_exponent=0.503-0.11*log(1.0-fppii)

      rh=2.16*(real(npep)**(v_exponent))
     &    +0.26*real(net_charge)
     &    -0.29*(real(npep)**(0.5))
      v_flory=log(rh/2.16)/log(real(npep))

      turn=0.0

      turn=num_ala*a_prop+num_cys*c_prop
     &        +num_asp*d_prop+num_glu*e_prop
     &        +num_phe*f_prop+num_gly*g_prop
     &        +num_his*h_prop+num_ile*i_prop
     &        +num_lys*k_prop+num_leu*l_prop
     &        +num_met*m_prop+num_asn*n_prop
     &        +num_pro*p_prop+num_gln*q_prop
     &        +num_arg*r_prop+num_ser*s_prop
     &        +num_thr*t_prop+num_val*v_prop
     &        +num_trp*w_prop+num_tyr*y_prop
      turn=turn/real(npep)

      sum_fb_hps=0.0

      sum_fb_hps=num_ala*fb_hps_a+num_cys*fb_hps_c
     &        +num_asp*fb_hps_d+num_glu*fb_hps_e
     &        +num_phe*fb_hps_f+num_gly*fb_hps_g
     &        +num_his*fb_hps_h+num_ile*fb_hps_i
     &        +num_lys*fb_hps_k+num_leu*fb_hps_l
     &        +num_met*fb_hps_m+num_asn*fb_hps_n
     &        +num_pro*fb_hps_p+num_gln*fb_hps_q
     &        +num_arg*fb_hps_r+num_ser*fb_hps_s
     &        +num_thr*fb_hps_t+num_val*fb_hps_v
     &        +num_trp*fb_hps_w+num_tyr*fb_hps_y
      sum_fb_hps=sum_fb_hps/real(npep)

      fb_hps(num_set)=sum_fb_hps
      scaling_exp(num_set)=v_flory
      AA_prop(num_set)=turn

      goto 10
11    close(2)


c Read input sequences for second set of proteins, using fasta format.

      open (2,file='PS_sequences_224.fasta',status='old',err=70)
      num_set2=0

20    continue
      read(2,22,err=12,end=400) code_inp
      length=len(code_inp)
      if(code_inp(1:1).eq.'>') goto 20

c  convert any lower case letter to upcase.
    
      do i=1,length
         num=ichar(code_inp(i:i))
         if (num.ge.97.and.num.le.122) code_inp(i:i) = char(num-32)
      enddo

c  Determine sequence length.

      j=0
      do i=1,length
      if (code_inp(i:i).eq.' ') goto 3
         j=j+1
         code(j)=code_inp(i:i)
3     continue
      enddo 
      npep=j

300   read(2,22,err=12,end=400) code_inp
      length=len(code_inp)
      if(code_inp(1:1).eq.'>') goto 400
      if(npep.ge.10000) goto 300

c  Add to sequence length for multiline entries.
    
      do i=1,length
         num=ichar(code_inp(i:i))
         if (num.ge.97.and.num.le.122) code_inp(i:i) = char(num-32)
      enddo

      j=npep
      do i=1,length
      if (code_inp(i:i).eq.' ') goto 4
         j=j+1
         code(j)=code_inp(i:i)
4     continue
      enddo
      npep=j
      
      goto 300

400   continue
      if(npep.gt.10000) goto 20
      num_set2=num_set2+1

      fb_hps(num_set+num_set2)=0.0
      scaling_exp(num_set+num_set2)=0.0
      AA_prop(num_set+num_set2)=0.0

c calculate v and turn propensity

      num_ala=0
      num_cys=0
      num_asp=0
      num_glu=0
      num_phe=0
      num_gly=0
      num_his=0
      num_ile=0
      num_lys=0
      num_leu=0
      num_met=0
      num_asn=0
      num_pro=0
      num_gln=0
      num_arg=0
      num_ser=0
      num_thr=0
      num_val=0
      num_trp=0
      num_tyr=0

      DO J=1,NPEP
         IF (CODE(J).EQ.'A') THEN
            num_ala=num_ala+1
         endif
         IF (CODE(J).EQ.'C') THEN
            num_cys=num_cys+1
         endif
         IF (CODE(J).EQ.'D') THEN
            num_asp=num_asp+1
         endif
         IF (CODE(J).EQ.'E') THEN
            num_glu=num_glu+1
         endif
         IF (CODE(J).EQ.'F') THEN
            num_phe=num_phe+1
         endif
         IF (CODE(J).EQ.'G') THEN
            num_gly=num_gly+1
         endif
         IF (CODE(J).EQ.'H') THEN
            num_his=num_his+1
         endif
         IF (CODE(J).EQ.'I') THEN
            num_ile=num_ile+1
         endif
         IF (CODE(J).EQ.'K') THEN
            num_lys=num_lys+1
         endif
         IF (CODE(J).EQ.'L') THEN
            num_leu=num_leu+1
         endif
         IF (CODE(J).EQ.'M') THEN
            num_met=num_met+1
         endif
         IF (CODE(J).EQ.'N') THEN
            num_asn=num_asn+1
         endif
         IF (CODE(J).EQ.'P') THEN
            num_pro=num_pro+1
         endif
         IF (CODE(J).EQ.'Q') THEN
            num_gln=num_gln+1
         endif
         IF (CODE(J).EQ.'R') THEN
            num_arg=num_arg+1
         endif
         IF (CODE(J).EQ.'S') THEN
            num_ser=num_ser+1
         endif
         IF (CODE(J).EQ.'T') THEN
            num_thr=num_thr+1
         endif
         IF (CODE(J).EQ.'V') THEN
            num_val=num_val+1
         endif
         IF (CODE(J).EQ.'W') THEN
            num_trp=num_trp+1
         endif
         IF (CODE(J).EQ.'Y') THEN
            num_tyr=num_tyr+1
         endif
      enddo

      net_charge=abs(num_asp+num_glu
     &                  -num_lys-num_arg)

      sum_ppii=0.0

      sum_ppii=num_ala*pppiia_h+num_cys*pppiic_h
     &        +num_asp*pppiid_h+num_glu*pppiie_h
     &        +num_phe*pppiif_h+num_gly*pppiig_h
     &        +num_his*pppiih_h+num_ile*pppiii_h
     &        +num_lys*pppiik_h+num_leu*pppiil_h
     &        +num_met*pppiim_h+num_asn*pppiin_h
     &        +num_pro*pppiip_h+num_gln*pppiiq_h
     &        +num_arg*pppiir_h+num_ser*pppiis_h
     &        +num_thr*pppiit_h+num_val*pppiiv_h
     &        +num_trp*pppiiw_h+num_tyr*pppiiy_h

      fppii=sum_ppii/real(npep)
      v_exponent=0.503-0.11*log(1.0-fppii)

      rh=2.16*(real(npep)**(v_exponent))
     &    +0.26*real(net_charge)
     &    -0.29*(real(npep)**(0.5))
      v_flory=log(rh/2.16)/log(real(npep))

      turn=0.0

      turn=num_ala*a_prop+num_cys*c_prop
     &        +num_asp*d_prop+num_glu*e_prop
     &        +num_phe*f_prop+num_gly*g_prop
     &        +num_his*h_prop+num_ile*i_prop
     &        +num_lys*k_prop+num_leu*l_prop
     &        +num_met*m_prop+num_asn*n_prop
     &        +num_pro*p_prop+num_gln*q_prop
     &        +num_arg*r_prop+num_ser*s_prop
     &        +num_thr*t_prop+num_val*v_prop
     &        +num_trp*w_prop+num_tyr*y_prop
      turn=turn/real(npep)

      sum_fb_hps=0.0

      sum_fb_hps=num_ala*fb_hps_a+num_cys*fb_hps_c
     &        +num_asp*fb_hps_d+num_glu*fb_hps_e
     &        +num_phe*fb_hps_f+num_gly*fb_hps_g
     &        +num_his*fb_hps_h+num_ile*fb_hps_i
     &        +num_lys*fb_hps_k+num_leu*fb_hps_l
     &        +num_met*fb_hps_m+num_asn*fb_hps_n
     &        +num_pro*fb_hps_p+num_gln*fb_hps_q
     &        +num_arg*fb_hps_r+num_ser*fb_hps_s
     &        +num_thr*fb_hps_t+num_val*fb_hps_v
     &        +num_trp*fb_hps_w+num_tyr*fb_hps_y
      sum_fb_hps=sum_fb_hps/real(npep)

      fb_hps(num_set+num_set2)=sum_fb_hps
      scaling_exp(num_set+num_set2)=v_flory
      AA_prop(num_set+num_set2)=turn

      goto 20
12    close(2)

      call UTEST(AA_prop,rank,num_set,num_set2,U,Z)
c      write(*,*)'The U-value is ',U
c      write(*,*)'The z-score is ',Z

      z_score=-z/sqrt(2.0)
      erf=0.0
      erfc=0.0
      
      call efunc(z_score,erf,erfc)

c output = num_AA_props, one-tail p-value
c      write(*,*)'one-tail p-value = ',0.5*erfc

      p_value=0.5*erfc

      if(p_value.lt.min_p_value) then
         min_p_value=p_value
         best_AA_prop=num_AA_props
      endif

      write(*,*)num_AA_props,p_value

      endif

      goto 80

800   continue
15    close(8)


      num_AA_props=num_AA_props+1
      call UTEST(scaling_exp,rank,num_set,num_set2,U,Z)
      z_score=-z/sqrt(2.0)
      erf=0.0
      erfc=0.0
      call efunc(z_score,erf,erfc)
      p_value=0.5*erfc
      if(p_value.lt.min_p_value) then
         min_p_value=p_value
         best_AA_prop=num_AA_props
      endif
      write(*,*)num_AA_props,p_value

      num_AA_props=num_AA_props+1
      call UTEST(fb_hps,rank,num_set,num_set2,U,Z)
      z_score=-z/sqrt(2.0)
      erf=0.0
      erfc=0.0
      call efunc(z_score,erf,erfc)
      p_value=0.5*erfc
      if(p_value.lt.min_p_value) then
         min_p_value=p_value
         best_AA_prop=num_AA_props
      endif
      write(*,*)num_AA_props,p_value

      write(*,*)'Best AA prop = ',best_AA_prop
      write(*,*)'min p-value = ',min_p_value

      end

C     ..................................................................
C
C        SUBROUTINE UTEST
C
C        PURPOSE
C           TEST WHETHER TWO INDEPENDENT GROUPS ARE FROM THE SAME
C           POPULATION BY MEANS OF MANN-WHITNEY U-TEST
C
C        USAGE
C           CALL UTEST(A,R,N1,N2,U,Z,IER)
C
C        DESCRIPTION OF PARAMETERS
C           A  - INPUT VECTOR OF CASES CONSISTING OF TWO INDEPENDENT
C                GROUPS. SMALLER GROUP PRECEDES LARGER GROUP. LENGTH
C                IS N1+N2.
C           R  - OUTPUT VECTOR OF RANKS. SMALLEST VALUE IS RANKED 1,
C                LARGEST IS RANKED N. TIES ARE ASSIGNED AVERAGE OF TIED
C                RANKS. LENGTH IS N1+N2.
C           N1 - NUMBER OF CASES IN SMALLER GROUP
C           N2 - NUMBER OF CASES IN LARGER GROUP
C           U  - STATISTIC USED TO TEST HOMOGENEITY OF THE TWO
C                GROUPS (OUTPUT)
C           Z  - MEASURE OF SIGNIFICANCE OF U IN TERMS OF NORMAL
C                DISTRIBUTION (OUTPUT)
C           IER- 0, IF NO ERROR.
C              - 1, IF ALL VALUES OF ONE GROUP ARE TIED.
C
C        REMARKS
C           Z IS SET TO ZERO IF N2 IS LESS THAN 20
C
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
C           RANK
C           TIE
C
C        METHOD
C           DESCRIBED IN S. SIEGEL, 'NONPARAMETRIC STATISTICS FOR THE
C           BEHAVIORAL SCIENCES', MCGRAW-HILL, NEW YORK, 1956,
C           CHAPTER 6
C
C     ..................................................................
C
      SUBROUTINE UTEST(A,R,N1,N2,U,Z)
      DIMENSION A(1000),R(1000)
C
C        RANK SCORES FROM BOTH GROUP TOGETHER IN ASCENDING ORDER, AND
C        ASSIGN TIED OBSERVATIONS AVERAGE OF TIED RANKS
C
      N=N1+N2
      CALL RANK(A,R,N)
      Z=0.0
C
C        SUM RANKS IN LARGER GROUP
C
      R2=0.0
      NP=N1+1

c fixed for deleted f77 feature
c      DO 10 I=NP,N
c   10 R2=R2+R(I)

      DO I=NP,N
        R2=R2+R(I)
      ENDDO

C
C        CALCULATE U
C
      FNX=N1*N2
      FN=N
      FN2=N2
      UP=FNX+FN2*((FN2+1.0)/2.0)-R2
      U=FNX-UP

c fixed for deleted f77 feature
c      IF(UP-U) 20,30,30
c   20 U=UP

      if ((UP-U).lt.(0)) U=UP

      if ((UP-U).ge.(0)) then

C        TEST FOR N2 LESS THAN 20

c fixed for deleted f77 feature
c   30 IF(N2-20) 80,40,40

      if ((N2-20).lt.0) goto 80

      if ((N2-20).ge.0) KT=1

C
C        COMPUTE STANDARD DEVIATION
C
c   40 KT=1

      endif

      CALL TIE(R,N,KT,TS)

c fixed for deleted f77 feature
c      IF(TS) 50,60,50

      if (TS.eq.(0.0)) S=SQRT(FNX*(FN+1.0)/12.0)

c      if (TS.lt.(0.0).or.TS.gt.(0.0)) then

c fixed for deleted f77 feature
c   50 IF (TS-(FN*FN*FN-FN)/12)52,51,52

      if ((TS-(FN*FN*FN-FN)/12).eq.(0.0)) then
         write(*,*)'error - ALL VALUES OF ONE GROUP ARE TIED'
         stop
      else
         S=SQRT((FNX/(FN*(FN-1.0)))*(((FN*FN*FN-FN)/12.0)-TS))
      endif

c   51 IER=1
c      GO TO 80
c   52 S=SQRT((FNX/(FN*(FN-1.0)))*(((FN*FN*FN-FN)/12.0)-TS))
c      GO TO 70
c   60 S=SQRT(FNX*(FN+1.0)/12.0)


C
C        COMPUTE Z
C
c   70 Z=(U-FNX*0.5)/S
c   80 RETURN

      Z=(U-FNX*0.5)/S
80    RETURN

      END

C     ..................................................................
C
C        SUBROUTINE RANK
C
C        PURPOSE
C           RANK A VECTOR OF VALUES
C
C        USAGE
C           CALL RANK(A,R,N)
C
C        DESCRIPTION OF PARAMETERS
C           A - INPUT VECTOR OF N VALUES
C           R - OUTPUT VECTOR OF LENGTH N. SMALLEST VALUE IS RANKED 1,
C               LARGEST IS RANKED N. TIES ARE ASSIGNED AVERAGE OF TIED
C               RANKS
C           N - NUMBER OF VALUES
C
C        REMARKS
C           NONE
C
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
C           NONE
C
C        METHOD
C           VECTOR IS SEARCHED FOR SUCCESSIVELY LARGER ELEMENTS. IF TIES
C           OCCUR, THEY ARE LOCATED AND THEIR RANK VALUE COMPUTED.
C           FOR EXAMPLE, IF 2 VALUES ARE TIED FOR SIXTH RANK, THEY ARE
C           ASSIGNED A RANK OF 6.5 (=(6+7)/2)
C
C     ..................................................................
C
      SUBROUTINE RANK(A,R,N)
      DIMENSION A(1000),R(1000)
C
C        INITIALIZATION
C
c fixed for deleted f77 feature
c      DO 10 I=1,N
c   10 R(I)=0.0

      DO I=1,N
         R(I)=0.0
      ENDDO

C
C        FIND RANK OF DATA
C
      DO 100 I=1,N
C
C        TEST WHETHER DATA POINT IS ALREADY RANKED
C

c fixed for deleted f77 feature
c      IF(R(I)) 20, 20, 100
C
C        DATA POINT TO BE RANKED
C

      if (R(I).gt.(0.0)) goto 100

   20 SMALL=0.0
      EQUAL=0.0
      X=A(I)
      DO 50 J=1,N

c fixed for deleted f77 feature
c      IF(A(J)-X) 30, 40, 50
C        COUNT NUMBER OF DATA POINTS WHICH ARE SMALLER
C
C

      if ((A(J)-X).lt.(0.0)) then
         SMALL=SMALL+1.0
         goto 50
      endif

c   30 SMALL=SMALL+1.0
c      GO TO 50
C
C        COUNT NUMBER OF DATA POINTS WHICH ARE EQUAL
C
c   40 EQUAL=EQUAL+1.0
c      R(J)=-1.0

      if ((A(J)-X).eq.(0.0)) then
         EQUAL=EQUAL+1.0
         R(J)=-1.0
      endif

      if ((A(J)-X).gt.(0.0)) goto 50

   50 CONTINUE

C
C        TEST FOR TIE
C

c fixed for deleted f77 feature
c      IF(EQUAL-1.0) 60, 60, 70
C
C        STORE RANK OF DATA POINT WHERE NO TIE
C
c   60 R(I)=SMALL+1.0
c      GO TO 100

      if ((EQUAL-1.0).le.(0.0)) then
         R(I)=SMALL+1.0
         goto 100
      endif

      if ((EQUAL-1.0).gt.(0.0)) then
C
C        CALCULATE RANK OF TIED DATA POINTS
C
   70 P=SMALL + (EQUAL + 1.0)*0.5
      DO 90 J=I,N

c fixed for deleted f77 feature
c      IF(R(J)+1.0) 90, 80, 90
c   80 R(J)=P

      if ((R(J)+1.0).eq.(0.0)) R(J)=P
   90 CONTINUE
      endif

  100 CONTINUE
      RETURN
      END

C     ..................................................................
C
C        SUBROUTINE TIE
C
C        PURPOSE
C           CALCULATE CORRECTION FACTOR DUE TO TIES
C
C        USAGE
C           CALL TIE(R,N,KT,T)
C
C        DESCRIPTION OF PARAMETERS
C           R  - INPUT VECTOR OF RANKS OF LENGTH N CONTAINING VALUES
C                1 TO N
C           N  - NUMBER OF RANKED VALUES
C           KT - INPUT CODE FOR CALCULATION OF CORRECTION FACTOR
C                      1   SOLVE EQUATION 1
C                      2   SOLVE EQUATION 2
C           T  - CORRECTION FACTOR (OUTPUT)
C                    EQUATION 1   T=SUM(CT**3-CT)/12
C                    EQUATION 2   T=SUM(CT*(CT-1)/2)
C                  WHERE CT IS THE NUMBER OF OBSERVATIONS TIED FOR A
C                        GIVEN RANK
C
C        REMARKS
C           NONE
C
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
C           NONE
C
C        METHOD
C           VECTOR IS SEARCHED FOR SUCCESSIVELY LARGER RANKS. TIES ARE
C           COUNTED AND CORRECTION FACTOR 1 OR 2 SUMMED.
C
C     ..................................................................
C
      SUBROUTINE TIE(R,N,KT,T)
      DIMENSION R(1000)
C
C        INITIALIZATION
C
      T=0.0
      Y=0.0
    5 X=1.0E38
      IND=0
C
C        FIND NEXT LARGEST RANK
C
      DO 30 I=1,N

c fixed for deleted f77 feature
c      IF(R(I)-Y) 30,30,10
c   10 IF(R(I)-X) 20,30,30
c   20 X=R(I)
c      IND=IND+1

      if ((R(I)-Y).gt.(0.0)) then
      if ((R(I)-X).lt.(0.0)) then
         X=R(I)
         IND=IND+1
      endif
      endif

   30 CONTINUE

C
C        IF ALL RANKS HAVE BEEN TESTED, RETURN
C
c      IF(IND) 90,90,40
c   40 Y=X
c      CT=0.0

      if (IND.le.(0)) goto 90
      if (IND.gt.(0)) Y=X
      CT=0.0

C
C        COUNT TIES
C
      DO 60 I=1,N

c fixed for deleted f77 feature
c      IF(R(I)-X) 60,50,60
c   50 CT=CT+1.0

      if ((R(I)-X).eq.(0.0)) CT=CT+1.0

   60 CONTINUE
C
C        CALCULATE CORRECTION FACTOR
C
c      IF(CT) 70,5,70
c   70 IF(KT-1) 75,80,75
c   75 T=T+CT*(CT-1.)/2.0
c      GO TO 5
c   80 T=T+(CT*CT*CT-CT)/12.0
c      GO TO 5

      if (CT.eq.(0.0)) goto 5

      if ((KT-1).eq.(0)) then
         T=T+(CT*CT*CT-CT)/12.0
         goto 5
      endif
      T=T+CT*(CT-1.)/2.0
      goto 5


   90 RETURN
      END

c
c From:
c An accurate solution FORTRAN algorithm for the erf and related error functions
c B.D. Brewster Geelen
c Advances in Engineering Software 18 (1993) 67-71
c
c CALC. OF ERROR FUNCTION ERF, ERFC FROM SERIES AND ASYMPTOTIC EXPANSIONS.
c
      SUBROUTINE EFUNC(X,ERF,ERFC)
      DOUBLE PRECISION X,E,PI,PKON,AN,FACT,A1,A2,A3,A4,A5,
     & ERF,ERFC,X2,F12,BT1,BT2,BT3,BT4,BT5,BT6,BT7,BT8,BT9,
     & SR,B1,B2,B3,B4

      E=2.718281828459045D0
      PI=4.0D0*DATAN(1.0D0)
      PKON=2.D0/(DSQRT(PI))
c
c SELECT EQUATION FOR VALUE OF X.
c

      IF(X.LT.(3.6)) GOTO 100
      GOTO 400
c
C *********************************************************
c
c      SUMMATION SOLUTION FOR X LT. 3.6
c

100   N=60 
      AN=-1.0D0
      FACT=1.0D0
      A5=0.0D0

      DO 200, J=1,N
      AN=AN+1.0D0
      FACT=FACT*AN
      IF(AN.EQ.0) FACT=1.0D0
      A1=(-1.0D0)**AN
      A2=(2.0D0*AN)+1.0D0
      A3=A1*X**A2
      A4=FACT*((2.0D0*AN)+1.0D0)
      A5=(A3/A4)+A5

      ERF=PKON*A5
      ERFC=1.0D0-ERF
200   CONTINUE
      GOTO 600

C *********************************************************
C
c       SERIES SOLUTION FOR X GT. 3.6
c

400   X2=2.0D0*X
      F12=132.0D0*3628800.0D0
      BT1=2.0D0/(X2**2.0D0)
      BT2=24.0D0/(2.0D0*(X2**4.0D0))
      BT3=720.0D0/(6.0D0*(X2**6.0D0))
      BT4=40320.0D0/(24.0D0*(X2**8.0D0))
      BT5=3628800.0D0/(120.0D0*(X2**10.0D0))
      BT6=F12/(720.0D0*(X2**12.0D0))
      BT7=(F12*182.0D0)/(5040.0D0*(X2**14.0D0))
      BT8=(F12*240.0D0*182.0D0)/(40320.0D0*(X2**16.0D0)) 
      BT9=(F12*306.0D0*240.0D0*182.0D0)/(362880.0D0*(X2**18.0D0))

      SR=1.0D0-BT1+BT2-BT3+BT4-BT5+BT6-BT7+BT8-BT9

      B1=X*X 
      B2=E**(-B1)
      B3=X*(DSQRT(PI))
      B4=B2/B3

      ERF=1.0D0-B4*SR
      ERFC=B4*SR

600   RETURN
      END