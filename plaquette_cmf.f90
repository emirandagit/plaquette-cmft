! Four-spin plaquette cluster mean-field solver for J=5/2 tetragonal layers.
! Energies: meV; temperature: K; magnetic field: T.
module plaquette_cmf
  use, intrinsic :: iso_fortran_env, only: real64,int64
  implicit none
  private
  integer, parameter :: dp=real64, ns=4, local_dim=6, dim=local_dim**ns
  integer, parameter :: max_nc=80, max_periods=8, max_links=16
  integer, parameter :: krylov=160, restarts=8
  real(dp), parameter :: pi=acos(-1.0_dp), kb=0.08617333262_dp
  real(dp), parameter :: mu_b=0.05788381806_dp, lande=6.0_dp/7.0_dp
  integer, parameter :: site_sub(ns)=[1,2,2,1]
  integer, parameter :: bond_i(4)=[1,1,2,3], bond_j(4)=[2,3,4,4]

  type, public :: parameters
    real(dp) :: j0=1.0_dp, j1=0.25_dp, j2=0.2022542486_dp
    real(dp) :: k0=0.0_dp, k1=0.0_dp, k2=0.0_dp
    real(dp) :: delta=0.8_dp, b20=0.0_dp, b40=0.0_dp, b44=0.0_dp
    real(dp) :: temperature=0.0_dp, theta_deg=0.0_dp
    real(dp) :: b_start=0.0_dp, b_stop=5.0_dp, b_step=0.25_dp
    real(dp) :: mixing=0.2_dp, tolerance=1.0e-7_dp
    integer :: n_periods=3, periods(max_periods)=[3,4,10,1,1,1,1,1]
    integer :: n_starts=10, max_iter=400
    integer :: start_ids(10)=[1,2,3,4,5,6,7,8,9,10]
    integer :: n_quad_starts=2, quad_start_ids(5)=[1,2,3,4,5]
    integer :: random_starts=0, random_seed=1729, seed_stride=1
    integer :: keep_branches=6
    integer :: anderson_depth=0
    logical :: dual_mixing=.false.
    real(dp) :: random_scale=0.45_dp
    logical :: reverse_sweep=.false.
    character(len=8) :: inplane='x'
    character(len=128) :: output_prefix='cmf'
  end type

  type, public :: solution
    real(dp) :: f=huge(1.0_dp), m(3,2,max_nc)=0.0_dp
    real(dp) :: qm(2,2,max_nc)=0.0_dp
    real(dp) :: residual=huge(1.0_dp)
    integer :: nc=0, seed=0, qseed=0, variant=0, method=0, iterations=0
    logical :: converged=.false.
  end type

  integer :: digit(ns,dim), p6(ns), nlinks(dim), link_row(max_links,dim)
  real(dp) :: diag0(dim), link_value(max_links,dim)
  real(dp) :: magnetic_m(0:5), up(0:5), down(0:5)
  real(dp) :: quad_up(0:5), quad_down(0:5)
  real(dp) :: coupl_j1, coupl_j2, coupl_k1, coupl_k2
  real(dp) :: rng_scale
  integer :: rng_seed

  public :: initialize, solve_branch, field_vector, order_wavevector
  public :: quadrupole_wavevector, net_moment, noncollinearity
  public :: dim, max_nc, max_periods, dp, pi, lande

contains

  subroutine add_link(col,row,value)
    integer, intent(in) :: col,row
    real(dp), intent(in) :: value
    integer :: k
    k=nlinks(col)+1
    if (k>max_links) error stop 'Increase max_links'
    nlinks(col)=k
    link_row(k,col)=row
    link_value(k,col)=value
  end subroutine

  subroutine initialize(par)
    type(parameters), intent(in) :: par
    integer :: col,s,d,i,j,a,b,t,row
    real(dp) :: z,jj,cf,coef
    jj=2.5_dp*(2.5_dp+1.0_dp)
    coupl_j1=par%j1; coupl_j2=par%j2
    coupl_k1=par%k1; coupl_k2=par%k2
    rng_seed=par%random_seed
    rng_scale=par%random_scale
    p6=[1,6,36,216]
    do d=0,5
      z=real(d,dp)-2.5_dp
      magnetic_m(d)=z
      up(d)=0.0_dp; down(d)=0.0_dp
      if (d<5) up(d)=sqrt(jj-z*(z+1.0_dp))
      if (d>0) down(d)=sqrt(jj-z*(z-1.0_dp))
      ! Qx+iQy={J+,Jz}/2, hence <m+1|Q+|m>=(m+1/2)<m+1|J+|m>.
      quad_up(d)=(z+0.5_dp)*up(d)
      quad_down(d)=(z-0.5_dp)*down(d)
    end do
    nlinks=0; diag0=0.0_dp
    do col=1,dim
      do s=1,ns
        digit(s,col)=mod((col-1)/p6(s),6)
      end do
      do s=1,ns
        d=digit(s,col); z=magnetic_m(d)
        cf=par%b20*(3.0_dp*z*z-jj) &
          +par%b40*(35.0_dp*z**4-(30.0_dp*jj-25.0_dp)*z*z &
          +3.0_dp*jj*jj-6.0_dp*jj)
        diag0(col)=diag0(col)+cf
        if (d<=1) then
          coef=1.0_dp
          do t=d,d+3
            coef=coef*up(t)
          end do
          row=col+4*p6(s)
          call add_link(col,row,0.5_dp*par%b44*coef)
        end if
        if (d>=4) then
          coef=1.0_dp
          do t=d,d-3,-1
            coef=coef*down(t)
          end do
          row=col-4*p6(s)
          call add_link(col,row,0.5_dp*par%b44*coef)
        end if
      end do
      do b=1,4
        i=bond_i(b); j=bond_j(b)
        a=digit(i,col); d=digit(j,col)
        diag0(col)=diag0(col)+par%j0*par%delta*magnetic_m(a)*magnetic_m(d)
        if (a<5 .and. d>0) then
          row=col+p6(i)-p6(j)
          call add_link(col,row,0.5_dp*(par%j0*up(a)*down(d) &
            +par%k0*quad_up(a)*quad_down(d)))
        end if
        if (a>0 .and. d<5) then
          row=col-p6(i)+p6(j)
          call add_link(col,row,0.5_dp*(par%j0*down(a)*up(d) &
            +par%k0*quad_down(a)*quad_up(d)))
        end if
      end do
    end do
  end subroutine

  subroutine apply_h(field,quad_field,vec,out)
    real(dp), intent(in) :: field(3,2),quad_field(2,2)
    complex(dp), intent(in) :: vec(dim)
    complex(dp), intent(out) :: out(dim)
    integer :: col,s,d,k,row,a
    real(dp) :: diagonal
    complex(dp) :: amp
    out=cmplx(0.0_dp,0.0_dp,dp)
    do col=1,dim
      diagonal=diag0(col)
      do s=1,ns
        a=site_sub(s); d=digit(s,col)
        diagonal=diagonal+field(3,a)*magnetic_m(d)
      end do
      amp=vec(col)
      out(col)=out(col)+diagonal*amp
      do k=1,nlinks(col)
        row=link_row(k,col)
        out(row)=out(row)+link_value(k,col)*amp
      end do
      do s=1,ns
        a=site_sub(s); d=digit(s,col)
        if (d<5) then
          row=col+p6(s)
          out(row)=out(row)+0.5_dp*(up(d)* &
            cmplx(field(1,a),-field(2,a),dp)+quad_up(d)* &
            cmplx(quad_field(1,a),-quad_field(2,a),dp))*amp
        end if
        if (d>0) then
          row=col-p6(s)
          out(row)=out(row)+0.5_dp*(down(d)* &
            cmplx(field(1,a),field(2,a),dp)+quad_down(d)* &
            cmplx(quad_field(1,a),quad_field(2,a),dp))*amp
        end if
      end do
    end do
  end subroutine

  subroutine moments(vec,m,qm)
    complex(dp), intent(in) :: vec(dim)
    real(dp), intent(out) :: m(3,2)
    real(dp), intent(out) :: qm(2,2)
    integer :: col,s,a,d,row
    complex(dp) :: c
    m=0.0_dp; qm=0.0_dp
    do col=1,dim
      do s=1,ns
        a=site_sub(s); d=digit(s,col)
        m(3,a)=m(3,a)+magnetic_m(d)*abs(vec(col))**2
        if (d<5) then
          row=col+p6(s)
          c=conjg(vec(row))*vec(col)
          m(1,a)=m(1,a)+real(c,dp)*up(d)
          m(2,a)=m(2,a)+aimag(c)*up(d)
          qm(1,a)=qm(1,a)+real(c,dp)*quad_up(d)
          qm(2,a)=qm(2,a)+aimag(c)*quad_up(d)
        end if
      end do
    end do
    m=m/2.0_dp
    qm=qm/2.0_dp
  end subroutine

  subroutine ground_cluster(field,quad_field,psi,e,m,qm,ok)
    real(dp), intent(in) :: field(3,2),quad_field(2,2)
    complex(dp), intent(inout) :: psi(dim)
    real(dp), intent(out) :: e,m(3,2),qm(2,2)
    logical, intent(out) :: ok
    complex(dp), allocatable :: q(:,:),w(:)
    real(dp) :: alpha(krylov),beta(krylov),dd(krylov),ee(krylov)
    real(dp), allocatable :: z(:,:)
    real(dp) :: work(2*krylov),res,bprev,normv
    complex(dp) :: overlap
    integer :: cycle,k,j,pass,info
    external :: dstev
    allocate(q(dim,krylov),w(dim),z(krylov,krylov))
    ok=.false.; normv=sqrt(real(dot_product(psi,psi),dp))
    if (normv<1.0e-20_dp) then
      do j=1,dim
        psi(j)=cmplx(cos(0.017_dp*j),sin(0.023_dp*j),dp)
      end do
      normv=sqrt(real(dot_product(psi,psi),dp))
    end if
    psi=psi/normv
    ! Keep overlap with other symmetry sectors when levels cross as B changes.
    do j=1,dim
      psi(j)=psi(j)+0.05_dp/sqrt(real(dim,dp))* &
        cmplx(cos(0.0137_dp*j*j),sin(0.0191_dp*j*j),dp)
    end do
    psi=psi/sqrt(real(dot_product(psi,psi),dp))
    do cycle=1,restarts
      q(:,1)=psi; bprev=0.0_dp
      do k=1,krylov
        call apply_h(field,quad_field,q(:,k),w)
        if (k>1) w=w-bprev*q(:,k-1)
        alpha(k)=real(dot_product(q(:,k),w),dp)
        w=w-alpha(k)*q(:,k)
        do pass=1,2
          do j=1,k
            overlap=dot_product(q(:,j),w)
            w=w-overlap*q(:,j)
          end do
        end do
        beta(k)=sqrt(max(0.0_dp,real(dot_product(w,w),dp)))
        if (mod(k,8)==0 .or. beta(k)<1.0e-13_dp .or. k==krylov) then
          dd(1:k)=alpha(1:k)
          if (k>1) ee(1:k-1)=beta(1:k-1)
          call dstev('V',k,dd,ee,z,krylov,work,info)
          if (info/=0) error stop 'LAPACK DSTEV failed'
          res=beta(k)*abs(z(k,1))
          if (res<1.0e-10_dp .or. beta(k)<1.0e-13_dp .or. k==krylov) then
            psi=cmplx(0.0_dp,0.0_dp,dp)
            do j=1,k
              psi=psi+z(j,1)*q(:,j)
            end do
            normv=sqrt(real(dot_product(psi,psi),dp))
            psi=psi/normv
            e=dd(1)
            if (res<1.0e-10_dp .or. beta(k)<1.0e-13_dp) then
              ok=.true.
              call moments(psi,m,qm)
              return
            end if
          end if
        end if
        if (k==krylov) exit
        q(:,k+1)=w/beta(k)
        bprev=beta(k)
      end do
    end do
    ! Dense exact diagonalization is a robust fallback for a difficult spectrum.
    call dense_cluster(field,quad_field,0.0_dp,psi,e,m,qm)
    ok=.true.
  end subroutine

  subroutine dense_cluster(field,quad_field,temperature,psi,free_energy,m,qm)
    real(dp), intent(in) :: field(3,2),quad_field(2,2),temperature
    complex(dp), intent(inout) :: psi(dim)
    real(dp), intent(out) :: free_energy,m(3,2),qm(2,2)
    complex(dp), allocatable :: h(:,:),v(:,:),work(:),basis(:),column(:)
    real(dp), allocatable :: eig(:),rwork(:)
    integer, allocatable :: iwork(:),isuppz(:)
    real(dp) :: tiny_abstol,weights(dim),partition,mp(3,2),qp(2,2),thermal_energy
    integer :: col,i,nfound,info,ihi
    external :: zheevr
    allocate(h(dim,dim),v(dim,dim),work(64*dim),basis(dim),column(dim))
    allocate(eig(dim),rwork(64*dim),iwork(16*dim),isuppz(2*dim))
    basis=cmplx(0.0_dp,0.0_dp,dp)
    do col=1,dim
      basis(col)=cmplx(1.0_dp,0.0_dp,dp)
      call apply_h(field,quad_field,basis,column)
      h(:,col)=column
      basis(col)=cmplx(0.0_dp,0.0_dp,dp)
    end do
    ihi=dim
    if (temperature<=0.0_dp) ihi=1
    tiny_abstol=0.0_dp
    call zheevr('V','I','U',dim,h,dim,0.0_dp,0.0_dp,1,ihi, &
      tiny_abstol,nfound,eig,v,dim,isuppz,work,size(work), &
      rwork,size(rwork),iwork,size(iwork),info)
    if (info/=0 .or. nfound/=ihi) error stop 'LAPACK ZHEEVR failed'
    psi=v(:,1)
    if (temperature<=0.0_dp) then
      free_energy=eig(1)
      call moments(psi,m,qm)
    else
      thermal_energy=kb*temperature
      weights(1:ihi)=exp(-(eig(1:ihi)-eig(1))/thermal_energy)
      partition=sum(weights(1:ihi))
      free_energy=eig(1)-thermal_energy*log(partition)
      m=0.0_dp; qm=0.0_dp
      do i=1,ihi
        if (weights(i)/partition<1.0e-14_dp) cycle
        call moments(v(:,i),mp,qp)
        m=m+weights(i)*mp/partition
        qm=qm+weights(i)*qp/partition
      end do
    end if
  end subroutine

  subroutine mean_field(par,nc,m,qm,h,hq)
    type(parameters), intent(in) :: par
    integer, intent(in) :: nc
    real(dp), intent(in) :: m(3,2,nc)
    real(dp), intent(in) :: qm(2,2,nc)
    real(dp), intent(out) :: h(3,2,nc),hq(2,2,nc)
    integer :: p,a,other,p1,p2,p3,p4
    do p=1,nc
      p1=modulo(p-2,nc)+1; p2=modulo(p,nc)+1
      p3=modulo(p-3,nc)+1; p4=modulo(p+1,nc)+1
      do a=1,2
        other=3-a
        h(:,a,p)=2.0_dp*par%j0*m(:,other,p)
        h(3,a,p)=2.0_dp*par%j0*par%delta*m(3,other,p)
        h(:,a,p)=h(:,a,p)+par%j1*(m(:,a,p1)+m(:,a,p2)) &
          +par%j2*(m(:,a,p3)+m(:,a,p4))
        hq(:,a,p)=2.0_dp*par%k0*qm(:,other,p) &
          +par%k1*(qm(:,a,p1)+qm(:,a,p2)) &
          +par%k2*(qm(:,a,p3)+qm(:,a,p4))
      end do
    end do
  end subroutine

  subroutine field_vector(par,b,bvec)
    type(parameters), intent(in) :: par
    real(dp), intent(in) :: b
    real(dp), intent(out) :: bvec(3)
    real(dp) :: theta,xy
    theta=par%theta_deg*pi/180.0_dp
    xy=b*sin(theta)
    if (trim(par%inplane)=='diag') then
      bvec=[xy/sqrt(2.0_dp),xy/sqrt(2.0_dp),b*cos(theta)]
    else
      bvec=[xy,0.0_dp,b*cos(theta)]
    end if
  end subroutine

  real(dp) function pseudo_random(variant,nc,p,component)
    integer, intent(in) :: variant,nc,p,component
    integer(int64) :: x
    x=modulo(int(rng_seed,int64)+104729_int64*variant+ &
      1009_int64*nc+7919_int64*p+65537_int64*component,2147483647_int64)
    x=modulo(48271_int64*x,2147483647_int64)
    pseudo_random=2.0_dp*real(x,dp)/2147483647.0_dp-1.0_dp
  end function

  subroutine rotate_xy(v,angle)
    real(dp), intent(inout) :: v(2)
    real(dp), intent(in) :: angle
    real(dp) :: old(2)
    old=v
    v(1)=cos(angle)*old(1)-sin(angle)*old(2)
    v(2)=sin(angle)*old(1)+cos(angle)*old(2)
  end subroutine

  subroutine seed_order(nc,seed,qseed,variant,bvec,m,qm)
    integer, intent(in) :: nc,seed,qseed,variant
    real(dp), intent(in) :: bvec(3)
    real(dp), intent(out) :: m(3,2,nc)
    real(dp), intent(out) :: qm(2,2,nc)
    real(dp) :: axis(3),e1(3),e2(3),l(3),phase,q,qk,base,phi,phiq
    real(dp) :: uniform_seed
    real(dp) :: jitter(3),qjitter(2),length_m
    integer :: p,a,signum
    axis=[0.0_dp,0.0_dp,1.0_dp]
    if (sqrt(sum(bvec*bvec))>1.0e-12_dp) axis=bvec/sqrt(sum(bvec*bvec))
    uniform_seed=min(1.8_dp,0.08_dp+ &
      0.05_dp*max(0.0_dp,sqrt(sum(bvec*bvec))-5.0_dp))
    e2=[-axis(2),axis(1),0.0_dp]
    if (sqrt(sum(e2*e2))<1.0e-10_dp) e2=[0.0_dp,1.0_dp,0.0_dp]
    e2=e2/sqrt(sum(e2*e2))
    e1=[e2(2)*axis(3)-e2(3)*axis(2), &
        e2(3)*axis(1)-e2(1)*axis(3), &
        e2(1)*axis(2)-e2(2)*axis(1)]
    if (coupl_j2<=coupl_j1/4.0_dp) then
      q=pi
    else
      q=acos(-coupl_j1/(4.0_dp*coupl_j2))
    end if
    q=2.0_dp*pi*real(nint(q*nc/(2.0_dp*pi)),dp)/real(nc,dp)
    qk=0.0_dp
    if (coupl_k1>0.0_dp) qk=pi
    if (coupl_k2>0.0_dp) then
      if (abs(coupl_k1/(4.0_dp*coupl_k2))<=1.0_dp) &
        qk=acos(-coupl_k1/(4.0_dp*coupl_k2))
    end if
    qk=2.0_dp*pi*real(nint(qk*nc/(2.0_dp*pi)),dp)/real(nc,dp)
    do p=1,nc
      phase=q*real(p-1,dp)
      select case(seed)
      case(1,2)
        signum=1
        if (seed==2) signum=-1
        l=1.8_dp*(e1*cos(phase)+real(signum,dp)*e2*sin(phase))
      case(3,4)
        signum=1
        if (modulo(p-1,3)==2) signum=-1
        l=1.8_dp*real(signum,dp)*e1
        if (seed==4) l=1.8_dp*real(signum,dp)*[0.0_dp,0.0_dp,1.0_dp]
      case(5,6)
        signum=1
        if (modulo(p-1,4)>=2) signum=-1
        l=1.8_dp*real(signum,dp)*e1
        if (seed==6) l=1.8_dp*real(signum,dp)*[0.0_dp,0.0_dp,1.0_dp]
      case(7)
        l=1.8_dp*e1
        if (modulo(p,2)==0) l=-l
      case(8)
        l=1.8_dp*e1
      case(9)
        l=1.8_dp*e2
      case(10)
        l=0.0_dp
      case default
        base=real(seed-10,dp)
        l=1.8_dp*(e1*cos(2.0_dp*pi*base*(p-1)/nc) &
          +e2*sin(2.0_dp*pi*base*(p-1)/nc))
      end select
      if (seed==10) then
        m(:,1,p)=1.8_dp*axis
        m(:,2,p)=1.8_dp*axis
      else
        m(:,1,p)=l+uniform_seed*axis
        m(:,2,p)=-l+uniform_seed*axis
      end if
      do a=1,2
        length_m=sqrt(sum(m(:,a,p)**2))
        if (length_m>2.4_dp) m(:,a,p)=m(:,a,p)*2.4_dp/length_m
      end do
      select case(qseed)
      case(1)
        do a=1,2
          qm(:,a,p)=0.4_dp*m(1:2,a,p)*m(3,a,p)
        end do
      case(2)
        qm(:,1,p)=1.5_dp*[cos(qk*(p-1)),sin(qk*(p-1))]
        qm(:,2,p)=qm(:,1,p)
      case(3)
        qm(:,1,p)=m(1:2,1,p)
        qm(:,2,p)=qm(:,1,p)
      case(4)
        qm(:,1,p)=[1.5_dp,0.0_dp]
        qm(:,2,p)=qm(:,1,p)
      case(5)
        qm(:,:,p)=0.0_dp
      case default
        error stop 'Invalid quadrupole seed'
      end select
    end do
    if (variant>0) then
      phi=pi*pseudo_random(variant,nc,0,1)
      phiq=pi*pseudo_random(variant,nc,0,2)
      do p=1,nc
        do a=1,2
          call rotate_xy(m(1:2,a,p),phi)
          call rotate_xy(qm(:,a,p),phiq)
        end do
        do a=1,3
          jitter(a)=1.5_dp*pseudo_random(variant,nc,p,2+a)
        end do
        do a=1,2
          qjitter(a)=1.5_dp*pseudo_random(variant,nc,p,5+a)
        end do
        m(:,1,p)=m(:,1,p)+rng_scale*jitter
        m(:,2,p)=m(:,2,p)-rng_scale*jitter
        qm(:,1,p)=qm(:,1,p)+rng_scale*qjitter
        qm(:,2,p)=qm(:,2,p)+rng_scale*qjitter
        do a=1,2
          length_m=sqrt(sum(m(:,a,p)**2))
          if (length_m>2.4_dp) m(:,a,p)=m(:,a,p)*2.4_dp/length_m
        end do
      end do
    end if
  end subroutine

  subroutine accelerated_step(x,fx,fhist,rhist,nhist,beta,xnext)
    real(dp), intent(in) :: x(:),fx(:),fhist(:,:),rhist(:,:),beta
    integer, intent(in) :: nhist
    real(dp), intent(out) :: xnext(size(x))
    real(dp) :: matrix(9,9),rhs(9,1),scale,regularizer
    real(dp) :: combined(size(x))
    integer :: pivot(9),i,j,info
    external :: dgesv
    xnext=x+beta*(fx-x)
    if (nhist<2) return
    matrix=0.0_dp; rhs=0.0_dp
    do i=1,nhist
      do j=1,nhist
        matrix(i,j)=dot_product(rhist(:,i),rhist(:,j))/size(x)
      end do
    end do
    scale=maxval(abs(matrix(1:nhist,1:nhist)))
    regularizer=1.0e-8_dp*scale+1.0e-16_dp
    do i=1,nhist
      matrix(i,i)=matrix(i,i)+regularizer
      matrix(i,nhist+1)=1.0_dp
      matrix(nhist+1,i)=1.0_dp
    end do
    rhs(nhist+1,1)=1.0_dp
    call dgesv(nhist+1,1,matrix,9,pivot,rhs,9,info)
    if (info/=0) return
    if (maxval(abs(rhs(1:nhist,1)))>10.0_dp) return
    combined=0.0_dp
    do i=1,nhist
      combined=combined+rhs(i,1)*fhist(:,i)
    end do
    if (maxval(abs(combined))>10.0_dp) return
    xnext=(1.0_dp-beta)*x+beta*combined
  end subroutine

  subroutine solve_branch(par,nc,seed,qseed,variant,bvec,initial_m,initial_q,result)
    type(parameters), intent(in) :: par
    integer, intent(in) :: nc,seed,qseed,variant
    real(dp), intent(in) :: bvec(3)
    real(dp), intent(in), optional :: initial_m(3,2,nc),initial_q(2,2,nc)
    type(solution), intent(out) :: result
    real(dp) :: m(3,2,nc),newm(3,2,nc),h(3,2,nc),effective(3,2)
    real(dp) :: qm(2,2,nc),newq(2,2,nc),hq(2,2,nc),effectiveq(2,2)
    real(dp) :: cluster_f(nc),residual,local_m(3,2),local_q(2,2),correction
    real(dp), allocatable :: x(:),fx(:),xnext(:),fhist(:,:),rhist(:,:)
    complex(dp), allocatable :: psi(:,:)
    integer :: p,a,it,j,nhist,nvar
    logical :: ok
    result=solution()
    result%nc=nc; result%seed=seed; result%qseed=qseed
    result%variant=variant
    if (par%anderson_depth>0) result%method=1
    residual=huge(1.0_dp)
    if (present(initial_m) .and. present(initial_q)) then
      m=initial_m; qm=initial_q
    else
      call seed_order(nc,seed,qseed,variant,bvec,m,qm)
    end if
    allocate(psi(dim,nc))
    nhist=0; nvar=10*nc
    if (par%anderson_depth>0) &
      allocate(x(nvar),fx(nvar),xnext(nvar), &
        fhist(nvar,par%anderson_depth),rhist(nvar,par%anderson_depth))
    do p=1,nc
      do j=1,dim
        psi(j,p)=cmplx(cos(0.017_dp*j+0.31_dp*p+seed), &
          sin(0.023_dp*j+0.19_dp*p+seed),dp)
      end do
    end do
    do it=1,par%max_iter
      call mean_field(par,nc,m,qm,h,hq)
      !$omp parallel do default(shared) &
      !$omp& private(p,a,effective,effectiveq,local_m,local_q,ok) schedule(dynamic)
      do p=1,nc
        effective=h(:,:,p)
        effectiveq=hq(:,:,p)
        do a=1,2
          effective(:,a)=effective(:,a)-lande*mu_b*bvec
        end do
        if (par%temperature<=0.0_dp) then
          call ground_cluster(effective,effectiveq,psi(:,p),cluster_f(p), &
            local_m,local_q,ok)
          if (.not.ok) error stop 'Ground-state diagonalization failed'
        else
          call dense_cluster(effective,effectiveq,par%temperature,psi(:,p), &
            cluster_f(p),local_m,local_q)
        end if
        newm(:,:,p)=local_m
        newq(:,:,p)=local_q
      end do
      !$omp end parallel do
      residual=max(maxval(abs(newm-m)),maxval(abs(newq-qm)))
      if (residual<par%tolerance) then
        correction=0.0_dp
        do p=1,nc
          do a=1,2
            correction=correction+dot_product(h(:,a,p),newm(:,a,p))
            correction=correction+dot_product(hq(:,a,p),newq(:,a,p))
          end do
        end do
        result%f=(sum(cluster_f)-correction)/(4.0_dp*nc)
        result%m(:,:,1:nc)=newm
        result%qm(:,:,1:nc)=newq
        result%residual=residual
        result%iterations=it
        result%converged=.true.
        return
      end if
      if (par%anderson_depth>0) then
        x(1:6*nc)=reshape(m,[6*nc])
        x(6*nc+1:nvar)=reshape(qm,[4*nc])
        fx(1:6*nc)=reshape(newm,[6*nc])
        fx(6*nc+1:nvar)=reshape(newq,[4*nc])
        if (nhist==par%anderson_depth) then
          fhist(:,1:nhist-1)=fhist(:,2:nhist)
          rhist(:,1:nhist-1)=rhist(:,2:nhist)
        else
          nhist=nhist+1
        end if
        fhist(:,nhist)=fx
        rhist(:,nhist)=fx-x
        call accelerated_step(x,fx,fhist,rhist,nhist,par%mixing,xnext)
        m=reshape(xnext(1:6*nc),shape(m))
        qm=reshape(xnext(6*nc+1:nvar),shape(qm))
      else
        m=m+par%mixing*(newm-m)
        qm=qm+par%mixing*(newq-qm)
      end if
    end do
    result%m(:,:,1:nc)=m
    result%qm(:,:,1:nc)=qm
    result%residual=residual
    result%iterations=par%max_iter
  end subroutine

  subroutine order_wavevector(result,q,amplitude)
    type(solution), intent(in) :: result
    real(dp), intent(out) :: q,amplitude
    complex(dp) :: fourier(3),phase
    real(dp) :: norm,best
    integer :: k,p,a,nc
    nc=result%nc
    q=0.0_dp; amplitude=0.0_dp; best=-1.0_dp
    do k=0,nc/2
      fourier=cmplx(0.0_dp,0.0_dp,dp)
      do p=1,nc
        phase=exp(cmplx(0.0_dp,-2.0_dp*pi*k*(p-1)/nc,dp))
        do a=1,3
          fourier(a)=fourier(a)+ &
            (result%m(a,1,p)-result%m(a,2,p))*phase/2.0_dp
        end do
      end do
      norm=sqrt(sum(abs(fourier)**2))/real(nc,dp)
      if (norm>best) then
        best=norm; q=real(k,dp)/real(nc,dp); amplitude=norm
      end if
    end do
  end subroutine

  subroutine quadrupole_wavevector(result,q,amplitude)
    type(solution), intent(in) :: result
    real(dp), intent(out) :: q,amplitude
    complex(dp) :: fourier(2),phase
    real(dp) :: norm,best
    integer :: k,p,a,nc
    nc=result%nc
    q=0.0_dp; amplitude=0.0_dp; best=-1.0_dp
    do k=0,nc/2
      fourier=cmplx(0.0_dp,0.0_dp,dp)
      do p=1,nc
        phase=exp(cmplx(0.0_dp,-2.0_dp*pi*k*(p-1)/nc,dp))
        do a=1,2
          fourier(a)=fourier(a)+ &
            (result%qm(a,1,p)+result%qm(a,2,p))*phase/2.0_dp
        end do
      end do
      norm=sqrt(sum(abs(fourier)**2))/real(nc,dp)
      if (norm>best) then
        best=norm; q=real(k,dp)/real(nc,dp); amplitude=norm
      end if
    end do
  end subroutine

  real(dp) function noncollinearity(result)
    type(solution), intent(in) :: result
    real(dp) :: u(3),v(3),cross(3),numerator,denominator
    integer :: p,pnext
    numerator=0.0_dp; denominator=0.0_dp
    do p=1,result%nc
      pnext=modulo(p,result%nc)+1
      u=(result%m(:,1,p)-result%m(:,2,p))/2.0_dp
      v=(result%m(:,1,pnext)-result%m(:,2,pnext))/2.0_dp
      cross=[u(2)*v(3)-u(3)*v(2),u(3)*v(1)-u(1)*v(3), &
        u(1)*v(2)-u(2)*v(1)]
      numerator=numerator+sqrt(sum(cross**2))
      denominator=denominator+sqrt(sum(u**2)*sum(v**2))
    end do
    noncollinearity=0.0_dp
    if (denominator>1.0e-12_dp) noncollinearity=numerator/denominator
  end function

  subroutine net_moment(result,m)
    type(solution), intent(in) :: result
    real(dp), intent(out) :: m(3)
    integer :: p
    m=0.0_dp
    do p=1,result%nc
      m=m+(result%m(:,1,p)+result%m(:,2,p))/(2.0_dp*result%nc)
    end do
  end subroutine

end module

program run_plaquette_cmf
  use plaquette_cmf
  implicit none
  integer, parameter :: bank_capacity=12
  type(parameters) :: par,par_linear
  type(solution) :: candidate,best_field
  type(solution) :: bank(bank_capacity,max_periods),next_bank(bank_capacity,max_periods)
  type(solution) :: endpoint_bank(bank_capacity,max_periods)
  type(solution), allocatable :: best_at(:)
  real(dp) :: b,bvec(3),q,amplitude,qq,qamp,net(3)
  integer :: bank_n(max_periods),next_n(max_periods)
  integer :: endpoint_n(max_periods)
  integer, allocatable :: attempts_at(:),converged_at(:)
  integer :: nfield,ib,ip,nc,seed,qseed,start_id,quad_id,random_start
  integer :: sweep,nsweeps,step,slot,unit1,unit2,unit3,unit4,p,a,ios
  integer :: variant
  logical :: fresh
  character(len=128) :: input_file
  character(len=16) :: subname(2)=[character(len=16) :: 'A','B']
  real(dp) :: j0,j1,j2,k0,k1,k2,delta,b20,b40,b44,temperature,theta_deg
  real(dp) :: b_start,b_stop,b_step,mixing,tolerance,random_scale
  integer :: n_periods,periods(max_periods),n_starts,max_iter,start_ids(10)
  integer :: n_quad_starts,quad_start_ids(5),random_starts,random_seed
  integer :: seed_stride,keep_branches,anderson_depth
  logical :: reverse_sweep,dual_mixing
  character(len=8) :: inplane
  character(len=128) :: output_prefix
  namelist /model/ j0,j1,j2,k0,k1,k2,delta,b20,b40,b44,temperature,theta_deg, &
    b_start,b_stop,b_step,mixing,tolerance,n_periods,periods,n_starts, &
    start_ids,n_quad_starts,quad_start_ids,max_iter,inplane,output_prefix, &
    random_starts,random_seed,random_scale,seed_stride,keep_branches, &
    anderson_depth,reverse_sweep,dual_mixing

  j0=par%j0; j1=par%j1; j2=par%j2; delta=par%delta
  k0=par%k0; k1=par%k1; k2=par%k2
  b20=par%b20; b40=par%b40; b44=par%b44
  temperature=par%temperature; theta_deg=par%theta_deg
  b_start=par%b_start; b_stop=par%b_stop; b_step=par%b_step
  mixing=par%mixing; tolerance=par%tolerance
  n_periods=par%n_periods; periods=par%periods
  n_starts=par%n_starts; start_ids=par%start_ids; max_iter=par%max_iter
  n_quad_starts=par%n_quad_starts; quad_start_ids=par%quad_start_ids
  random_starts=par%random_starts; random_seed=par%random_seed
  random_scale=par%random_scale; seed_stride=par%seed_stride
  keep_branches=par%keep_branches; reverse_sweep=par%reverse_sweep
  anderson_depth=par%anderson_depth
  dual_mixing=par%dual_mixing
  inplane=par%inplane; output_prefix=par%output_prefix
  input_file='cmf.nml'
  if (command_argument_count()>0) call get_command_argument(1,input_file)
  open(newunit=unit1,file=trim(input_file),status='old',action='read',iostat=ios)
  if (ios/=0) error stop 'Cannot open input namelist'
  read(unit1,nml=model,iostat=ios)
  close(unit1)
  if (ios/=0) error stop 'Cannot parse input namelist'
  par%j0=j0; par%j1=j1; par%j2=j2; par%delta=delta
  par%k0=k0; par%k1=k1; par%k2=k2
  par%b20=b20; par%b40=b40; par%b44=b44
  par%temperature=temperature; par%theta_deg=theta_deg
  par%b_start=b_start; par%b_stop=b_stop; par%b_step=b_step
  par%mixing=mixing; par%tolerance=tolerance
  par%n_periods=n_periods; par%periods=periods
  par%n_starts=n_starts; par%start_ids=start_ids; par%max_iter=max_iter
  par%n_quad_starts=n_quad_starts; par%quad_start_ids=quad_start_ids
  par%random_starts=random_starts; par%random_seed=random_seed
  par%random_scale=random_scale; par%seed_stride=seed_stride
  par%keep_branches=keep_branches; par%reverse_sweep=reverse_sweep
  par%anderson_depth=anderson_depth
  par%dual_mixing=dual_mixing
  par%inplane=inplane; par%output_prefix=output_prefix
  if (par%j0<=0.0_dp .or. par%j1<=0.0_dp .or. par%j2<=0.0_dp) &
    error stop 'All dipole exchange constants must be positive'
  if (par%delta<=0.0_dp .or. par%delta>=1.0_dp) &
    error stop 'Delta must lie strictly between 0 and 1'
  if (par%temperature<0.0_dp .or. par%b_step<=0.0_dp .or. &
    par%b_stop<par%b_start) error stop 'Invalid temperature or field grid'
  if (par%n_periods<1 .or. par%n_periods>max_periods) error stop 'Invalid n_periods'
  if (any(par%periods(1:par%n_periods)<1) .or. &
    any(par%periods(1:par%n_periods)>max_nc)) error stop 'Invalid period'
  if (par%n_starts<1 .or. par%n_starts>size(par%start_ids) .or. &
    par%max_iter<1 .or. par%mixing<=0.0_dp .or. &
    par%mixing>1.0_dp .or. par%tolerance<=0.0_dp) error stop 'Invalid solver settings'
  if (any(par%start_ids(1:par%n_starts)<1)) error stop 'Invalid start_ids'
  if (par%n_quad_starts<1 .or. par%n_quad_starts>size(par%quad_start_ids)) &
    error stop 'Invalid n_quad_starts'
  if (any(par%quad_start_ids(1:par%n_quad_starts)<1) .or. &
    any(par%quad_start_ids(1:par%n_quad_starts)>5)) &
    error stop 'Invalid quad_start_ids'
  if (par%random_starts<0 .or. par%random_scale<0.0_dp .or. &
    par%seed_stride<1 .or. par%keep_branches<1 .or. &
    par%keep_branches>bank_capacity) error stop 'Invalid search settings'
  if (par%anderson_depth<0 .or. par%anderson_depth>8) &
    error stop 'anderson_depth must be between 0 and 8'
  if (trim(par%inplane)/='x' .and. trim(par%inplane)/='diag') &
    error stop 'inplane must be x or diag'

  call initialize(par)
  par_linear=par
  par_linear%anderson_depth=0
  nfield=int((par%b_stop-par%b_start)/par%b_step+1.0e-9_dp)
  allocate(best_at(0:nfield),attempts_at(0:nfield),converged_at(0:nfield))
  best_at=solution(); attempts_at=0; converged_at=0
  open(newunit=unit3,file=trim(par%output_prefix)//'_branches.csv', &
    status='replace',action='write')
  write(unit3,'(a)') 'B_T,sweep,Nc,seed,qseed,variant,method,F_meV_per_spin,'// &
    'Qz_dip_rlu,dip_amplitude,Qz_quad_rlu,quad_amplitude,'// &
    'noncollinearity,residual,iterations'

  nsweeps=1
  if (par%reverse_sweep) nsweeps=2
  do sweep=1,nsweeps
    if (sweep==1) then
      bank_n=0
    else
      bank_n=endpoint_n
      bank=endpoint_bank
    end if
    do step=0,nfield
      ib=step
      if (sweep==2) ib=nfield-step
      b=par%b_start+par%b_step*ib
      call field_vector(par,b,bvec)
      best_field=solution()
      do ip=1,par%n_periods
        nc=par%periods(ip)
        next_n(ip)=0
        do slot=1,bank_n(ip)
          attempts_at(ib)=attempts_at(ib)+1
          if (par%dual_mixing) then
            call solve_branch(par_linear,nc,0,0,0,bvec, &
              bank(slot,ip)%m(:,:,1:nc),bank(slot,ip)%qm(:,:,1:nc),candidate)
          else
            call solve_branch(par,nc,0,0,0,bvec, &
              bank(slot,ip)%m(:,:,1:nc),bank(slot,ip)%qm(:,:,1:nc),candidate)
          end if
          call record_candidate(candidate)
        end do
        fresh=mod(step,par%seed_stride)==0 .or. step==nfield
        if (next_n(ip)==0) fresh=.true.
        if (fresh) call fresh_starts()
        if (next_n(ip)>0) then
          bank_n(ip)=next_n(ip)
          bank(1:bank_n(ip),ip)=next_bank(1:bank_n(ip),ip)
        end if
        write(*,'(a,i1,a,f7.2,a,i3,a,i3)') &
          'sweep=',sweep,' B=',b,' T Nc=',nc,' retained=',bank_n(ip)
      end do
      if (best_field%converged) then
        if (best_field%f<best_at(ib)%f) best_at(ib)=best_field
        call order_wavevector(best_field,q,amplitude)
        write(*,'(a,i1,a,f7.2,a,i3,a,f16.8,a,f6.3)') &
          'sweep=',sweep,' B=',b,' T Nc=',best_field%nc, &
          ' F=',best_field%f,' Qd=',q
      else
        write(*,'(a,i1,a,f7.2,a)') 'sweep=',sweep,' B=',b, &
          ' T: no converged branch'
      end if
      flush(unit3)
    end do
    if (sweep==1) then
      endpoint_n=bank_n
      endpoint_bank=bank
    end if
  end do
  close(unit3)

  open(newunit=unit1,file=trim(par%output_prefix)//'_summary.csv', &
    status='replace',action='write')
  open(newunit=unit2,file=trim(par%output_prefix)//'_moments.csv', &
    status='replace',action='write')
  open(newunit=unit4,file=trim(par%output_prefix)//'_quadrupoles.csv', &
    status='replace',action='write')
  write(unit1,'(a)') 'B_T,Nc,seed,qseed,variant,method,F_meV_per_spin,'// &
    'Qz_dip_rlu,dip_amplitude,Qz_quad_rlu,quad_amplitude,'// &
    'noncollinearity,Mx_muB,My_muB,Mz_muB,residual,iterations,'// &
    'converged_branches,attempted_branches'
  write(unit2,'(a)') 'B_T,Nc,plane,sublattice,Jx,Jy,Jz'
  write(unit4,'(a)') 'B_T,Nc,plane,sublattice,Qx,Qy'
  do ib=0,nfield
    if (.not.best_at(ib)%converged) error stop 'No branch converged at a field'
    b=par%b_start+par%b_step*ib
    call order_wavevector(best_at(ib),q,amplitude)
    call quadrupole_wavevector(best_at(ib),qq,qamp)
    call net_moment(best_at(ib),net)
    write(unit1,'(*(g0,:,","))') b,best_at(ib)%nc,best_at(ib)%seed, &
      best_at(ib)%qseed,best_at(ib)%variant,best_at(ib)%method,best_at(ib)%f, &
      q,amplitude,qq,qamp,noncollinearity(best_at(ib)), &
      lande*net(1),lande*net(2),lande*net(3), &
      best_at(ib)%residual,best_at(ib)%iterations, &
      converged_at(ib),attempts_at(ib)
    do p=1,best_at(ib)%nc
      do a=1,2
        write(unit2,'(*(g0,:,","))') b,best_at(ib)%nc,p,trim(subname(a)), &
          best_at(ib)%m(1,a,p),best_at(ib)%m(2,a,p),best_at(ib)%m(3,a,p)
        write(unit4,'(*(g0,:,","))') b,best_at(ib)%nc,p,trim(subname(a)), &
          best_at(ib)%qm(1,a,p),best_at(ib)%qm(2,a,p)
      end do
    end do
  end do
  close(unit1); close(unit2); close(unit4)

contains

  subroutine fresh_starts()
    integer :: base_index,quad_index
    do seed=1,par%n_starts
      start_id=par%start_ids(seed)
      do qseed=1,par%n_quad_starts
        quad_id=par%quad_start_ids(qseed)
        attempts_at(ib)=attempts_at(ib)+1
        call solve_branch(par,nc,start_id,quad_id,0,bvec,result=candidate)
        call record_candidate(candidate)
        if (par%dual_mixing .and. par%anderson_depth>0) then
          attempts_at(ib)=attempts_at(ib)+1
          call solve_branch(par_linear,nc,start_id,quad_id,0,bvec,result=candidate)
          call record_candidate(candidate)
        end if
      end do
    end do
    do random_start=1,par%random_starts
      base_index=1+mod(random_start-1,par%n_starts)
      quad_index=1+mod((random_start-1)/par%n_starts,par%n_quad_starts)
      start_id=par%start_ids(base_index)
      quad_id=par%quad_start_ids(quad_index)
      variant=random_start+1000*(sweep-1)+10000*ib
      attempts_at(ib)=attempts_at(ib)+1
      call solve_branch(par,nc,start_id,quad_id,variant,bvec,result=candidate)
      call record_candidate(candidate)
      if (par%dual_mixing .and. par%anderson_depth>0) then
        attempts_at(ib)=attempts_at(ib)+1
        call solve_branch(par_linear,nc,start_id,quad_id,variant,bvec,result=candidate)
        call record_candidate(candidate)
      end if
    end do
  end subroutine

  subroutine record_candidate(found)
    type(solution), intent(in) :: found
    real(dp) :: dipq,dipamp,quadq,quadamp
    if (.not.found%converged) return
    converged_at(ib)=converged_at(ib)+1
    if (found%f<best_field%f) best_field=found
    call insert_branch(found,next_bank(:,ip),next_n(ip))
    call order_wavevector(found,dipq,dipamp)
    call quadrupole_wavevector(found,quadq,quadamp)
    write(unit3,'(*(g0,:,","))') b,sweep,nc,found%seed,found%qseed, &
      found%variant,found%method,found%f,dipq,dipamp,quadq,quadamp, &
      noncollinearity(found),found%residual,found%iterations
    flush(unit3)
  end subroutine

  subroutine insert_branch(found,collection,nstored)
    type(solution), intent(in) :: found
    type(solution), intent(inout) :: collection(bank_capacity)
    integer, intent(inout) :: nstored
    integer :: i,worst
    do i=1,nstored
      if (same_branch(found,collection(i))) then
        if (found%f<collection(i)%f) collection(i)=found
        return
      end if
    end do
    if (nstored<par%keep_branches) then
      nstored=nstored+1
      collection(nstored)=found
      return
    end if
    worst=1
    do i=2,nstored
      if (collection(i)%f>collection(worst)%f) worst=i
    end do
    if (found%f<collection(worst)%f) collection(worst)=found
  end subroutine

  logical function same_branch(left,right)
    type(solution), intent(in) :: left,right
    real(dp) :: distance
    integer :: shift,swap,p,other,site
    same_branch=.false.
    if (left%nc/=right%nc) return
    if (abs(left%f-right%f)>1.0e-4_dp) return
    do shift=0,left%nc-1
      do swap=0,1
        distance=0.0_dp
        do p=1,left%nc
          other=modulo(p-1+shift,left%nc)+1
          do site=1,2
            distance=max(distance,maxval(abs(left%m(:,site,p)- &
              right%m(:,1+modulo(site-1+swap,2),other))))
            distance=max(distance,maxval(abs(left%qm(:,site,p)- &
              right%qm(:,1+modulo(site-1+swap,2),other))))
          end do
        end do
        if (distance<0.03_dp) then
          same_branch=.true.
          return
        end if
      end do
    end do
  end function

end program
