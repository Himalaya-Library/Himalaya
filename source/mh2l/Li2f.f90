!>    @brief C wrapper for complex dilogarithm
      subroutine li2c(re_in, im_in, re_out, im_out) bind(C, name="li2c_")
      implicit none
      double precision, intent(in) :: re_in, im_in
      double precision, intent(out) :: re_out, im_out
      double complex z, l, cdli2

      z = dcmplx(re_in, im_in)
      l = cdli2(z)

      re_out = dble(l)
      im_out = dimag(l)

      return
      end
