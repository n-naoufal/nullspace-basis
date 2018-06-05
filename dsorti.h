interface
    subroutine dsorti(which, apply, n, x1, x2)
        integer  n
        character*2  which
        logical  apply
        Double precision  x1(0:n-1)
        Double precision  x2(0:n-1)
    end subroutine dsorti
end interface


