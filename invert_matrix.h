/** 
 * \file invert_matrix.h
 * \brief Invert matrix through QR decomposition
 * \author Pascal Getreuer <getreuer@gmail.com>
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under, at your option, the terms of the GNU General Public License as 
 * published by the Free Software Foundation, either version 3 of the 
 * License, or (at your option) any later version, or the terms of the 
 * simplified BSD license.
 *
 * You should have received a copy of these licenses along this program. 
 * If not, see <http://www.gnu.org/licenses/> and
 * <http://www.opensource.org/licenses/bsd-license.html>.
 */

#ifndef _INVERT_MATRIX_H_
#define _INVERT_MATRIX_H_

int invert_matrix(double *inv_A, double *A, int N);

#endif /* _INVMAT_H_ */
