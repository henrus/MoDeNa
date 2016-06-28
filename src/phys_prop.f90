!auxilliary subroutines for calculation of state properties
module phys_prop
    implicit none
    private
    public volume_balance,dispress,min_film_thickness
contains
!********************************BEGINNING*************************************
!checks whether we are losing some mass or not
pure subroutine  volume_balance(y,vt,fs)
    use globals
    real(dp), intent(out) :: vt,fs
    real(dp), dimension(:), intent(in) :: y
    integer :: i,neq
    real(dp) :: vf,vs
    neq=size(y)
    vf=0
    vs=0
    do i=1,neq
        if (y(i)<strutFilmParameter*y(1)) then
            vf=vf+2*pi*dr*(0.5_dp+i-1)*y(i)*dr
        else
            vs=vs+2*pi*dr*(0.5_dp+i-1)*y(i)*dr
        endif
    enddo
    vs=vs+pi*y(neq)/3*(rd**2+rc*rd+rc**2)-pi*y(neq)*rd**2
    vt=vf+vs
    fs=vs/vt
end subroutine volume_balance
!***********************************END****************************************


!********************************BEGINNING*************************************
!disjoining pressure
pure subroutine dispress(h,dispr,dph)
    use globals
    real(dp), intent(in) :: h
    real(dp), intent(out) :: dispr !disjoining pressure
    real(dp), intent(out) :: dph !derivative of disjoining pressure
    dispr=bdp*((hdp/h)**(ndp-1)-(hdp/h)**(mdp-1))*(hdp-cdp)
    dph=(bdp*(hdp/h)**mdp*(hdp*mdp+cdp*(h-h*mdp))-&
        bdp*(hdp/h)**ndp*(hdp*ndp+cdp*(h-h*ndp)))/(h*hdp)
end subroutine dispress
!***********************************END****************************************


!********************************BEGINNING*************************************
!minimum film thickness and its distance from the center
pure subroutine min_film_thickness(y,hmin,hloc)
    use globals
    real(dp), dimension(:), intent(in) :: y
    real(dp), intent(out) :: hmin
    real(dp), intent(out) :: hloc
    hmin=minval(y,dim=1)
    hloc=minloc(y,dim=1)*dr
end subroutine min_film_thickness
!***********************************END****************************************
end module phys_prop
