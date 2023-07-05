program Model_Recursive
    implicit none
    integer :: i, j, n
    real :: new_x, new_y
    real, dimension (:,:), allocatable :: A, H, Y, TH, h_k, y_k !TH is parameter matrix (Theta) !Inputs
    real, dimension (:,:), allocatable :: tr_h_k, M1, TR_H, M2, M5, M6, M3, M4, M7, TH_N !Intermediates for New Theta
    real, dimension (:,:), allocatable :: E1, E2, E3, FE
    interface
        subroutine Matrix_Trans(In1, In2)
        implicit none
        integer :: nrow, ncol, i, j
        real, dimension (:,:), allocatable :: In1, In2
        end subroutine Matrix_Trans      
    
        subroutine Matrix_Multi(In1,In2,Output)
        implicit none
        integer :: Row1, Col, Col2, i, j, k
        real, dimension (:,:), allocatable :: In1, In2, Output
        end subroutine Matrix_Multi
    
        subroutine MatrixInverse(Input1, Output1)
        implicit none    
        integer :: nrow, ncol
        real, dimension(:,:), allocatable :: Input1, Output1
        integer :: i, j, k
        real :: X11, Xik
        end subroutine MatrixInverse   
        
        subroutine Matrix_Add(Input1, Input2, Output)
        implicit none
        integer :: i,j, nrow, ncol
        real, dimension(:,:), allocatable :: Input1, Input2, Output
        end subroutine Matrix_Add
        
        subroutine Matrix_Sub(Input1, Input2, Output)
        implicit none
        integer :: i,j, nrow, ncol
        real, dimension(:,:), allocatable :: Input1, Input2, Output
        end subroutine Matrix_Sub      
    end interface
    open(1, file = "E:\Fortran\Exercies\Recursive Model\XY.txt")
    open(2,file = "E:\Fortran\Exercies\Recursive Model\test.txt")
    open(3, file = "E:\Fortran\Exercies\Recursive Model\THETA.txt")
    n = 0
    do while (.true.)
        read(1, *, iostat = i) 
        if (i /= 0) exit  ! Exit the loop when end of file is reached
        n = n + 1
    end do
    print *, n
    close(1)
    open(1, file = "E:\Fortran\Exercies\Recursive Model\XY.txt")
    allocate (A(n+1,2), H(n,2), Y(n,1), h_k(1,2), y_k(1,1), TH(2,1))
    do i = 1,n
        read(1,*) A(i,:)
    end do
    close(1)
    do i = 1,2
        read(3,*) TH(i,1)
    end do 
    close(2)
    close(3)
    do i = 1, n
        Y(i,1) = A(i,2)
    end do
    do i = 1, n
        H(i,1) = A(i,1)
        H(i,2) = 1
    end do 
    !Getting new data
    print *, "Enter new value of x"
    read (*,*) new_x
    print *, "Enter new value of y"
    read (*,*) new_y
    h_k(1,1) = new_x
    h_k(1,2) = 1
    y_k(1,1) = new_y
    
    !New Parameter Matrix
    call Matrix_Trans(h_k,tr_h_k)
    call Matrix_Multi(tr_h_k,h_k,M1)
    call Matrix_Trans(H, TR_H)
    call Matrix_Multi(TR_H, H, M2)
    call Matrix_Add(M1, M2, M5)
    call MatrixInverse(M5,M6)
    call Matrix_Multi(tr_h_k, y_k, M3)
    call Matrix_Multi(M2, TH, M4)
    call Matrix_Add(M3, M4, M7)
    call Matrix_Multi(M6, M7, TH_N)
    
    !Replacing the Parameter Text File
    open (3, file = "E:\Fortran\Exercies\Recursive Model\THETA.txt", status = 'replace')
    do i = 1,2
        write (3,*) TH_N(i,1)
    end do
    close(3)
    
    !Adding the new data to original
    A(n+1,1) = new_x
    A(n+1,2) = new_y
    open(1, file = "E:\Fortran\Exercies\Recursive Model\XY.txt", status = 'old', access = 'append')
    do j = 1,2
        write(1,'(2f5.2)', advance = 'no') A(n+1,j)
        write(1, '(a)', advance = 'no') char(9)
    end do
    close(1)    
    
    !Forecast Error
    call Matrix_Multi(h_k,TH,E1)
    call Matrix_Sub(y_k,E1,E2)
    call Matrix_Multi(tr_h_k, E2, E3)
    call Matrix_Multi(M6, E3, FE)
    do i = 1,2
        print *, FE(i,:)
    end do
    end program Model_Recursive
    
    
    subroutine Matrix_Trans(In1, In2)
    implicit none
    integer :: nrow, ncol, i, j
    real, dimension (:,:), allocatable :: In1, In2
    nrow = size(In1, dim = 1)
    ncol = size(In1, dim = 2)
    allocate (In2(ncol,nrow))

    do j = 1, ncol
        do i = 1,nrow
            In2(j,i) = In1(i,j)
        end do
    end do
    end subroutine Matrix_Trans

subroutine Matrix_Multi(In1,In2,Output)
    implicit none
    integer :: Row1, Col, Col2, i, j, k
    real, dimension (:,:), allocatable :: In1, In2, Output
    Row1 = size(In1, dim = 1)
    Col = size(In1, dim = 2)
    Col2 = size(In2, dim = 2)
    allocate (Output(Row1,Col2))

    do i = 1,Row1
        do j = 1,Col2
            Output(i,j) = 0.0
            do k = 1,Col
                Output(i,j) = Output(i,j) + In1(i,k)*In2(k,j)
            end do
        end do
    end do
    end subroutine Matrix_Multi
      
subroutine MatrixInverse(Input1, Output1)
    implicit none
    integer :: nrow, ncol
    real, dimension(:,:), allocatable :: Input1, Output1
    integer :: i, j, k
    real :: X11, Xik
    nrow = size(Input1, dim = 1)
    ncol = nrow
    allocate(Output1(nrow, ncol))
    
    !writing an Identity matrix
    do i = 1, nrow
        do j = 1, ncol
            if (i == j) then
                Output1(i,j) = 1
            else
                Output1(i,j) = 0
            end if
        end do
    end do
    
    do k = 1, nrow
        X11 = Input1(k,k)
        do j = 1, ncol
            Input1(k,j) = Input1(k,j)/X11
            Output1(k,j) = Output1(k,j)/X11
        end do
          
        do i = 1, nrow
            if (i /= k) then
                Xik = Input1(i,k)
                do j = 1, ncol
                  Input1(i,j) = Input1(i,j) - Xik * Input1(k,j)
                  Output1(i,j) = Output1(i,j) - Xik * Output1(k,j)
                end do
            end if
        end do
    end do
    end subroutine MatrixInverse 
    
    subroutine Matrix_Add(Input1, Input2, Output)
    implicit none
    integer :: i,j, nrow, ncol
    real, dimension(:,:), allocatable :: Input1, Input2, Output
    nrow = size(Input1, dim = 1)
    ncol = size(Input2, dim = 2)
    allocate (Output(nrow,ncol))
    do i = 1, nrow
        do j = 1, ncol
            Output(i,j) = Input1(i,j) + Input2(i,j)
        end do
    end do
    end subroutine Matrix_Add
    
    subroutine Matrix_Sub(Input1, Input2, Output)
    implicit none
    integer :: i,j, nrow, ncol
    real, dimension(:,:), allocatable :: Input1, Input2, Output
    nrow = size(Input1, dim = 1)
    ncol = size(Input2, dim = 2)
    allocate (Output(nrow,ncol))
    do i = 1, nrow
        do j = 1, ncol
            Output(i,j) = Input1(i,j) - Input2(i,j)
        end do
    end do
    end subroutine Matrix_Sub