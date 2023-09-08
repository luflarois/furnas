module modStatistic

   !! ## Module with some statistics functions
   !!
   !! ## BRAMS-FURNAS
   !!
   !! Author: Rodrigues, L.F.
   !!
   !! E-mail: <mailto:luizfrodrigues@protonmail.com>
   !!
   !! Date: 23March2023 15:27
   !!
   !! #####Version: 0.1.0
   !!
   !! ---
   !! **Full description**:
   !!
   !! Module with some statistics functions
   !!
   !! ** History**:
   !!
   !! --- 
   !! ** Licence **:
   !!
   !!  <img src="https://www.gnu.org/graphics/gplv3-127x51.png width="63">
   !!
   !!  This program is free software: you can redistribute it and/or modify
   !!  it under the terms of the GNU General Public License as published by
   !!  the  Free  Software  Foundation, either version 3 of the License, or
   !!  (at your option) any later version.
   !!
   !!  This program is distributed in the hope that it  will be useful, but
   !!  ** WITHOUT  ANY  WARRANTY **;  without  even  the   implied   warranty  of
   !!  **MERCHANTABILITY** or **FITNESS FOR A  PARTICULAR PURPOSE**.  See  the, GNU
   !!  GNU General Public License for more details.
   !!
   !!  You should have received a copy  of the GNU General  Public  License
   !!  along with this program.  If not, see [GNU Public License](https://www.gnu.org/licenses/gpl-3.0.html).
   !!
   !use, intrinsic :: iso_c_binding
   !use m
   implicit none
   !use modConstants, only: pi
   character(len=*), parameter :: p_source_name = 'modStatistic.F90' 
   !! Source code name 
   character(len=*), parameter :: p_module_name = 'modStatistic' 
   !! module name 
   real, parameter :: pi = 3.1415926
   real, parameter :: sqrt2pi = sqrt(2.0*pi)

   private
   public :: Mu, Sigma2, PFunction, quantil

contains
   !----------------------------------------------------------------------------------------------------------
   function Mu(values_in) result(average)
      !! ## Compute a average of N values
      !!
      !! Author: Rodrigues, L.F.
      !!
      !! E-mail: <mailto:luizfrodrigues@protonmail.com>
      !!
      !! Date: 23March2023 15:29
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! Compute a average of N values
      !! $\mu={\sum_{i=1}^{n}{{x_i}}\over{n}}$
      !!
      !! ** History**:
      !!
      !! --- 
      !! ** Licence **: Under the terms of the GNU General Public version 3
      !!   <img src="https://www.gnu.org/graphics/gplv3-127x51.png width="63">
      !!
   
      implicit none
      ! Parameters:
      character(len=*), parameter :: p_procedure_ame = 'Mu' 
      !! Function name
   
      ! Variables (input):
      real, intent(in) :: values_in(:)
   
      ! Local variables:
      real :: average
      !! function Output: average
      integer :: sizeOfValues, icnt

      sizeOfvalues = size(values_in)

      average = 0.0
      do icnt=1,sizeOfvalues
         average = average + values_in(icnt)
      enddo
      average = average/real(sizeOfvalues)

   end function Mu

   !---------------------------------------------------------------------------------------------------
   function Sigma2(average,values_in) result(variance)
      !! ## Compute a variance of N values
      !!
      !! Author: Rodrigues, L.F.
      !!
      !! E-mail: <mailto:luizfrodrigues@protonmail.com>
      !!
      !! Date: 23March2023 15:37
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! Compute a variance of N values
      !! $\sigma^{2}={{\sum_{i=1}^{n}(x_i-\mu)^2}\over{n}}$
      !!
      !! ** History**:
      !!
      !! --- 
      !! ** Licence **: Under the terms of the GNU General Public version 3
      !!   <img src="https://www.gnu.org/graphics/gplv3-127x51.png width="63">
      !!
   
      implicit none
      ! Parameters:
      character(len=*), parameter :: p_procedure_ame = 'Sigma2' 
   
      !! Function name
   
      ! Variables (input):
      real, intent(in) :: values_in(:)
      !! Input values
      real, intent(in) :: average
      !! Average of same values
   
      ! Local variables:
      real :: variance
      !! function Output: variance
      integer :: sizeOfValues, icnt
   
      sizeOfvalues = size(values_in)

      variance = 0.0
      do icnt=1,sizeOfvalues
         variance = variance + (values_in(icnt)-average)**2
      end do
      variance = variance/sizeOfvalues
   
   end function Sigma2

   function PFunction(percentil,average,variance) result(p_result)
      !! ## Compoute the P function of P
      !!
      !! Author: Rodrigues, L.F.
      !!
      !! E-mail: <mailto:luizfrodrigues@protonmail.com>
      !!
      !! Date: 23March2023 15:53
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! Compoute the P function of P
      !! $f(p,\mu,\sigma)={{1}\over{\sqrt{2\pi}\sigma}}e^{-{(p-\mu)2} \over{2\sigma^{2}}}$
      !!
      !! ** History**:
      !!
      !! --- 
      !! ** Licence **: Under the terms of the GNU General Public version 3
      !!   <img src="https://www.gnu.org/graphics/gplv3-127x51.png width="63">
      !!
   
      implicit none
      ! Parameters:
      character(len=*), parameter :: p_procedure_ame = 'PFunction' 
      !! Function name
   
      ! Variables (input):
      real, intent(in) :: percentil
      !! The percentil of interest P90, P50, etc (0.9, 0.5, etc)
      real, intent(in) :: average
      !! Average of values
      real, intent(in) :: variance
      !! Variance of values
   
      ! Local variables:
      real :: p_result
      !! function Output: p_result
      real :: p1,p2
   
      ! Code:
      p1 = 1/(sqrt2pi*sqrt(variance))
      p2 = (-(percentil-average)**2)/(2*variance)
      p_result = p1*exp(p2)
   
   end function PFunction

   function quantil(list,qnt) result(qtl)

      implicit none

      real(8), intent(in) :: list(:)
      ! Vector of values
      real, intent(in) :: qnt
      ! Quantile to return
      integer :: size_of_list, i
      real :: position, pint
      integer, allocatable :: index(:)
      !
      real :: qtl
      !Return the quantile

      size_of_list = size(list)
      allocate(index(size_of_list))

      position = qnt * size_of_list
      pint  = int(position)
      ! print *,position
      ! do i = 1, size_of_list
      !    print *, i, list(i)
      ! end do
      qtl = 0.0

      if(position - pint > 0.0) then 
         qtl = 0.5 * (list(pint)+list(pint+1))
      else 
         qtl = list(int(position))
      end if

   end function quantil

end module modStatistic