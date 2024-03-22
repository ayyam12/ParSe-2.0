c----------------------------------------------------------------
c
c Calculates from sequence the polymer scaling exponent, nu, and
c the beta turn propensity for the whole input sequence.
c
c STW 07/02/2020
c
c input file (i.e., proteome in fasta format) should be named sequences.fasta (line 172)
c     protein length range is defined on lines 243 and 246
c     sliding window size is defined on line 253
c input file of AA scales should be named aaindex1_new.txt (line 119)
c     number of scales to use is defined on line 69
c     identity of the scales are specified beginning on line 71
c output file of the scale-calculated properties of each 25-residue window is PCA_dataset.txt
c
c----------------------------------------------------------------

      program PCA_dataset_generation
      implicit none

      integer length,num,i,j,npep,jj,num_seq,k,k_scale(100),
     &        count,num_AA_props,k_num,window_size,
     &        num_ala,num_cys,num_asp,num_glu,num_phe,
     &        num_gly,num_his,num_ile,num_lys,num_leu,
     &        num_met,num_asn,num_pro,num_gln,num_arg,
     &        num_ser,num_thr,num_val,num_trp,num_tyr
      character code_inp*11000,code(11000)*3,title*4000,type(500)*2,
     &          remark*3,uniprotID(100)*3
      real pppiia_h,pppiip_h,pppiig_h,pppiic_h,pppiid_h,pppiie_h,
     &     pppiif_h,pppiih_h,pppiii_h,pppiik_h,pppiil_h,pppiim_h,
     &     pppiin_h,pppiiq_h,pppiir_h,pppiis_h,pppiit_h,pppiiv_h,
     &     pppiiw_h,pppiiy_h,slope,intercept,v_line,net_charge,
     &     sum_ppii,v_exponent,fppii,v_flory,rh,
     &     a_prop(600),c_prop(600),d_prop(600),
     &     e_prop(600),f_prop(600),g_prop(600),h_prop(600),
     &     i_prop(600),k_prop(600),l_prop(600),m_prop(600),
     &     n_prop(600),p_prop(600),q_prop(600),r_prop(600),
     &     s_prop(600),t_prop(600),v_prop(600),w_prop(600),
     &     y_prop(600),scale(600)


c intrinsic PPII bias measured in peptides by Hilser group
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


c specify which scales to test.

      k_num=33
c top 3 helix scales (discerning test vs null)
      k_scale(1)=186
      k_scale(2)=366
      k_scale(3)=442
c top 3 turn scales (discerning test vs null)
      k_scale(4)=62
      k_scale(5)=166
      k_scale(6)=228
c top 3 sheet scales (discerning test vs null)
      k_scale(7)=283
      k_scale(8)=225
      k_scale(9)=172
c top 3 aperiodic scales (discerning test vs null)
      k_scale(10)=105
      k_scale(11)=104
      k_scale(12)=107
c top 3 coil and loop scales (discerning test vs null)
      k_scale(13)=350
      k_scale(14)=188
      k_scale(15)=351
c top 3 flexibility scales (discerning test vs null)
      k_scale(16)=505
      k_scale(17)=144
      k_scale(18)=439
c top 3 size scales (discerning test vs null)
      k_scale(19)=28
      k_scale(20)=81
      k_scale(21)=112
c top 3 composition scales (discerning test vs null)
      k_scale(22)=469
      k_scale(23)=461
      k_scale(24)=462
c top 3 charge scales (discerning test vs null)
      k_scale(25)=89
      k_scale(26)=88
      k_scale(27)=146
c top 3 hydrophobic_structure scales (discerning test vs null)
      k_scale(28)=365
      k_scale(29)=400
      k_scale(30)=114
c top 3 hydrophobic_solution scales (discerning test vs null)
      k_scale(31)=434
      k_scale(32)=178
      k_scale(33)=179


c read the different amino acid properties in the amino acid index.

      num_AA_props=0
      open (8,file='aaindex1_new.txt',status='old')
80    continue
      read(8,*,err=15,end=800) remark
      if(remark.eq.'I  ') then
      num_AA_props=num_AA_props+1
      read(8,*)a_prop(num_AA_props),r_prop(num_AA_props),
     &         n_prop(num_AA_props),d_prop(num_AA_props),
     &         c_prop(num_AA_props),q_prop(num_AA_props),
     &         e_prop(num_AA_props),g_prop(num_AA_props),
     &         h_prop(num_AA_props),i_prop(num_AA_props)
      read(8,*)l_prop(num_AA_props),k_prop(num_AA_props),
     &         m_prop(num_AA_props),f_prop(num_AA_props),
     &         p_prop(num_AA_props),s_prop(num_AA_props),
     &         t_prop(num_AA_props),w_prop(num_AA_props),
     &         y_prop(num_AA_props),v_prop(num_AA_props)
      endif
      goto 80
800   continue
15    close(8)
c      write(*,*)num_AA_props
c      do jj=1,num_AA_props
c      write(*,*)a_prop(jj),r_prop(jj),
c     &         n_prop(jj),d_prop(jj),
c     &         c_prop(jj),q_prop(jj),
c     &         e_prop(jj),g_prop(jj),
c     &         h_prop(jj),i_prop(jj)
c      write(*,*)l_prop(jj),k_prop(jj),
c     &         m_prop(jj),f_prop(jj),
c     &         p_prop(jj),s_prop(jj),
c     &         t_prop(jj),w_prop(jj),
c     &         y_prop(jj),v_prop(jj)
c      stop
c      enddo

      count=0

      open (3,file='PCA_dataset.txt')
c      write(3,'("Seq#  Win#  nu",33i6)')(k_scale(j),j=1,k_num)

      write(3,'("Seq#  Win#  nu_model",
     &          "  helix.1  helix.2  helix.3",
     &          "  turn.1  turn.2  turn.3",
     &          "  sheet.1  sheet.2  sheet.2",
     &          "  aperiodic.1  aperiodic.2  aperiodic.3",
     &          "  coil.1  coil.2  coil.3",
     &          "  flexibility.1  flexibility.2  flexibility.3",
     &          "  size.1  size.2  size.3",
     &          "  composition.1  composition.2  composition.3",
     &          "  neg.charge  pos.charge  net.charge",
     &          "  hydropho_str.1  hydropho_str.2  hydropho_str.3",
     &          "  hydropho_sol.1  hydropho_sol.2  hydropho_sol.3")')

c Read in the sequences in fasta format

      open (2,file='sequences.fasta',status='old',err=70)
      goto 10
70    write (*,'("No input file present")')
      stop
      
22    format (A11000)
10    continue
      read(2,22,err=11,end=200) code_inp
      length=len(code_inp)

c save uniprot ID info on the primary sequence
      if(code_inp(1:1).eq.'>') then
      j=0
      do i=1,length
         j=j+1
         if (j.gt.75) goto 3
         uniprotID(j)=code_inp(i:i)
3     continue
      enddo

      goto 10
      endif

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
      if(code_inp(1:1).eq.'>') then
         backspace(2) 
         goto 200
      endif
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

c if a sequence is shorter than 25 residues, thus shorter than the sliding window, then ignore
      if(npep.lt.50)  goto 10

c if a sequence is longer than 10000, it may exceed the allocated memory for a primary sequence, then ignore
      if(npep.gt.500) goto 10

      num_seq=num_seq+1


c calculate v and other intrinsic properties

      window_size=25

      DO J=1,NPEP
         if (j.le.(NPEP-window_size+1)) then

      count=count+1
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

      DO JJ=J,J+window_size-1
         IF (CODE(JJ).EQ.'A') THEN
            num_ala=num_ala+1
         endif
         IF (CODE(JJ).EQ.'C') THEN
            num_cys=num_cys+1
         endif
         IF (CODE(JJ).EQ.'D') THEN
            num_asp=num_asp+1
         endif
         IF (CODE(JJ).EQ.'E') THEN
            num_glu=num_glu+1
         endif
         IF (CODE(JJ).EQ.'F') THEN
            num_phe=num_phe+1
         endif
         IF (CODE(JJ).EQ.'G') THEN
            num_gly=num_gly+1
         endif
         IF (CODE(JJ).EQ.'H') THEN
            num_his=num_his+1
         endif
         IF (CODE(JJ).EQ.'I') THEN
            num_ile=num_ile+1
         endif
         IF (CODE(JJ).EQ.'K') THEN
            num_lys=num_lys+1
         endif
         IF (CODE(JJ).EQ.'L') THEN
            num_leu=num_leu+1
         endif
         IF (CODE(JJ).EQ.'M') THEN
            num_met=num_met+1
         endif
         IF (CODE(JJ).EQ.'N') THEN
            num_asn=num_asn+1
         endif
         IF (CODE(JJ).EQ.'P') THEN
            num_pro=num_pro+1
         endif
         IF (CODE(JJ).EQ.'Q') THEN
            num_gln=num_gln+1
         endif
         IF (CODE(JJ).EQ.'R') THEN
            num_arg=num_arg+1
         endif
         IF (CODE(JJ).EQ.'S') THEN
            num_ser=num_ser+1
         endif
         IF (CODE(JJ).EQ.'T') THEN
            num_thr=num_thr+1
         endif
         IF (CODE(JJ).EQ.'V') THEN
            num_val=num_val+1
         endif
         IF (CODE(JJ).EQ.'W') THEN
            num_trp=num_trp+1
         endif
         IF (CODE(JJ).EQ.'Y') THEN
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

      fppii=sum_ppii/real(window_size)
      v_exponent=0.503-0.11*log(1.0-fppii)

c multiplier of 4 on the sequence is because window_size = 25
c we found that sequence calculated v_flory is length dependent for sequences shorter than ~100 residues
c at sufficient chain lengths (≥100 residues) v_flory depends on composition not length

      rh=2.16*(real(4*window_size)**(v_exponent))
     &    +0.26*real(4*net_charge)
     &    -0.29*(real(4*window_size)**(0.5))
      v_flory=log(rh/2.16)/log(real(4*window_size))

      do k=1,k_num
      do jj=1,num_AA_props
      if(jj.eq.k_scale(k)) then
      scale(k)=0.0

      scale(k)=num_ala*a_prop(jj)+num_cys*c_prop(jj)+
     &         num_asp*d_prop(jj)+num_glu*e_prop(jj)+
     &         num_phe*f_prop(jj)+num_gly*g_prop(jj)+
     &         num_his*h_prop(jj)+num_ile*i_prop(jj)+
     &         num_lys*k_prop(jj)+num_leu*l_prop(jj)+
     &         num_met*m_prop(jj)+num_asn*n_prop(jj)+
     &         num_pro*p_prop(jj)+num_gln*q_prop(jj)+
     &         num_arg*r_prop(jj)+num_ser*s_prop(jj)+
     &         num_thr*t_prop(jj)+num_val*v_prop(jj)+
     &         num_trp*w_prop(jj)+num_tyr*y_prop(jj)
      scale(k)=scale(k)/real(window_size)
      endif
      enddo
      enddo

      write(3,'(i10,i10,f10.4,33f10.4)')
     &  num_seq,count,v_flory,(scale(k),k=1,k_num)

      endif
      enddo

      goto 10
11    close(2)


      close(3)

      end