/*
	Copyright Michael McCrea, 2019
		M McCrea	mtm5@uw.edu
	This file is part of the MatrixArray quark for SuperCollider 3 and is free software:
	you can redistribute it and/or modify it under the terms of the GNU General
	Public License as published by the Free Software Foundation, either version 3
	of the License, or (at your option) any later version.
	This code is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
	without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
	See the GNU General Public License for more details.
	<http://www.gnu.org/licenses/>.

	A description of the functionality, along with attributions can be found in
	both the README and help documentation.
*/

MatrixArray {

	var <matrix, <rows, <cols, <isSquare;
	var flopped;

	*newClear { |rows = 1, cols = 1|
		^super.new.init(
			rows.collect{
				cols.collect{ nil }
			}
		);
	}

	*newIdentity { |size = 1|
		^super.new.init(
			size.collect{ |r|
				Array.fill(size, { |c|
					if (r == c) {1}{0};
				})
			}
		);
	}

	*newUnitriangular { |size = 1|
		^super.new.init(
			size.collect{ |r|
				Array.fill(size, { |c|
					if (r <= c) {1}{0};
				})
			}
		);
	}

	*fill { |rows, cols, func|
		^super.new.init(
			rows.collect{ |r|
				Array.fill(cols, { |c|
					func.(r, c)
				})
			}
		);
	}

	// returns MatrixArray from 2D array (array of rows)
	*with { |array|
		^super.new.init(array)
	}

	init { |array|
		if (array.rank != 2) {
			error(
				format(
					"[MatrixArray:-init] Array must be rank 2 to create a matrix. Received rank %",
					array.rank
				)
			).throw
		} {
			if (array.collect(_.size).every(_ == array[0].size).not) {
				error("[MatrixArray:-init] Row sizes are not identical").throw
			}
		};

		// initialize instance vars
		matrix = array;
		rows = array.size;
		cols = array[0].size;
		isSquare = rows == cols;
		flopped = nil;
	}

	// get the flopped matrix (transpose), computing it if needed
	flopped {
		^flopped ?? { flopped = matrix.flop }
	}
	transpose { ^this.flopped }

	// returns the elements/row/col objects from the matrix (not copies!)
	at { |row, col| ^matrix[row][col] }
	rowAt { |row| ^matrix[row] }
	colAt { |col| ^this.flopped[col] }

	// return the array that is the matrix (not a copy!)
	asArray { ^matrix }

	put { |row, col, val|
		matrix[row][col] = val;
		flopped = nil; // trigger re-calc of transpose on next request
	}

	// this is a destructive operation:
	// force values to zero that are within threshold distance (positive or negative)
	zeroWithin { |within = (-180.dbamp)|
		var item;

		rows.do{ |r|
			cols.do{ |c|
				item = matrix[r][c];
				matrix[r][c] = if(item.abs <= within, { 0 }, { item });
			}
		};

		flopped = nil; // trigger re-calc of transpose on next request
	}


	// returns an Array of multiplying with either a number or matrix
	// A(m,n) * B(n,r) = AB(m,r)
	* { |that|
		^if (that.isNumber) {
			this.mulNumber(that);
		} {
			this.mulMatrix(that);
		};
	}

	// returns an Array of multiplying with a matrix
	// A(m,n) * B(n,r) = AB(m,r)
	// mArray can be a MatrixArray or (2D) Array
	mulMatrix { |mArray|
		var mulArr, mulArrFlopped, result;

		mulArr = mArray.asArray; // allow aMatrix as MatrixArray or Array

		if (cols == mulArr.size) { // this.cols == mArray.rows
			mulArrFlopped = mulArr.flop;
			result = Array.fill(rows, { Array.newClear(mulArr[0].size) });

			matrix.do { | row, i |
				mulArrFlopped.do { | tcol, j |
					result[i][j] = sum(row * tcol)
				}
			};

			^result
		} {
			format(
				"[MatrixArray:-mulMatrix] cols and rows don't fit. Matrix shapes: (%) (%)",
				matrix.shape, mulArr.shape
			).throw;
		};
	}

	// returns a new Array of multiplying with a number
	mulNumber { |aNumber|
		^matrix * aNumber;
	}

	// returns a new Array without row
	withoutRow { |row|
		var mtx = matrix.copy;
		mtx.removeAt(row);
		^mtx
	}
	// returns a new Array without col
	withoutCol { |col|
		var mtx = this.flopped.copy;
		mtx.removeAt(col);
		^mtx.flop
	}

	// convenience method: return an MatrixArray without row, col
	// returns MatrixArray for optimized use by cofactor
	withoutRowCol { |row, col|
		^MatrixArray.with(this.withoutRow(row)).removeCol(col);
	}

	// destructively removes row
	removeRow { |row|
		matrix.removeAt(row);
		this.init(matrix); // matrix changed, need to re-initialize state variables
	}
	// destructively removes col
	removeCol { |col|
		this.flopped.removeAt(col);
		matrix = flopped.flop;
		this.init(matrix); // matrix changed, need to re-initialize state variables
	}

	// return a new Array of arrays representing sub matrix within the matrix
	getSub { |rowStart = 0, colStart = 0, rowLength, colHeight|
		var w, h, sub, maxw, maxh;
		var rFrom, rTo;

		maxw = cols - rowStart;
		maxh = rows - colStart;

		w = rowLength ?? { maxw };
		h = colHeight ?? { maxh };

		if ((w > maxw) or: (h > maxh)) {
			format(
				"[MatrixArray:-getSub] Dimensions of requested sub-matrix exceed bounds: "
				"you asked for %x%, remaining space after starting index is %x%",
				rowLength, colHeight, maxw, maxh
			).throw
		};

		sub = Array.newClear(h);

		rFrom = rowStart;
		rTo = rFrom + rowLength - 1;

		(colStart .. (colStart + h - 1)).do{ |row, i|
			sub[i] = this.rowAt(row)[rFrom..rTo];
		};

		^sub
	}


	// return single value cofactor
	// NOTE: for efficiency, assumed to be square matrix (doesn't check)
	cofactor { |row, col|
		^((-1) ** (row + col)) * this.withoutRowCol(row, col).det
	}

	// return a new Array that is the gram of this matrix
	gram {
		// "cast" to a MatrixArray to force mulMatrix
		^MatrixArray.with(this.flopped) * this
	}

	// return a new Array that is the adjoint of this matrix
	adjoint {
		var adjoint = Array.fill(rows, { Array.newClear(cols) });

		rows.do{ |i|
			cols.do{ |j|
				adjoint[i][j] = this.cofactor(i, j);
			}
		};

		^adjoint.flop;
	}

	// return a new Array that is the inverse of this matrix
	inverse {
		var det = this.det;

		if (det != 0.0) {
			^this.adjoint / det;
		} {
			"[MatrixArray:-inverse] Matrix is singular.".throw
		};
	}

	// returns a new Array which is the pseudoInverse of the matrix
	pseudoInverse {
		var gram, grami, mul, muli;

		if (cols < rows) {
			gram = MatrixArray.with(this.gram);
			grami = MatrixArray.with(gram.inverse);
			^grami * this.flopped;
		} {
			mul = this * this.flopped;
			muli = MatrixArray.with(mul).inverse;
			^MatrixArray.with(this.flopped) * muli;
		};
	}

	// returns determinant via LU Decomposition
	// TODO: can determinant be calculated from Upper Triangular matrix only?
	// requires square matrix
	det {
		var i, p, k, j, ri, m, n;
		var im1, pm1, jm1, nm1;
		var det = 1.0;

		// TODO: investigate alternatives to deepCopy here: store temp matrix globally?
		m = matrix.deepCopy; // create a deep copy, operations happen in-place

		n = m.size;   // size of rows or cols (square)
		nm1 = n-1;
		ri = n.collect{ |i| i };

		// LU factorization.
		p = 1;
		while ({ p <= (nm1) }, {
			pm1 = p - 1;
			// Find pivot element.
			i = p + 1;
			while ({ i <= n }, {
				im1 = i - 1;
				if ( abs(m[ri[im1]][pm1]) > abs(m[ri[pm1]][pm1])) {
					var t;
					// Switch the index for the p-1 pivot row if necessary.
					t = ri[pm1];
					ri[pm1] = ri[im1];
					ri[im1] = t;
					det = det.neg;
				};
				i = i + 1;
			});

			if (m[ri[pm1]][pm1] == 0) { "[MatrixArray:-det] The matrix is singular.".throw };

			// Multiply the diagonal elements.
			det = det * m[ri[pm1]][pm1];

			// Form multiplier.
			i = p + 1;
			while ({ i <= n }, {
				im1 = i - 1;
				m[ri[im1]][pm1] = m[ri[im1]][pm1] / m[ri[pm1]][pm1];

				// Eliminate [p-1].
				j = p + 1;
				while( { j <= n }, {
					jm1 = j - 1;
					m[ri[im1]][jm1] = m[ri[im1]][jm1] - (m[ri[im1]][pm1] * m[ri[pm1]][jm1]);

					j = j + 1;
				});

				i = i + 1;
			});

			p = p + 1;
		});

		det = det * m[ri[nm1]][nm1];

		^det // return
	}

	// Mix coefficients with the matrix
	// (e.g. applying a tranform to ambisonic coefficients)
	// Returns a new Array of size this.rows
	mixCoeffs { |coeffs|

		if (coeffs.size != cols) {
			format(
				"[MatrixArray:-mixCoeffs] - coeffs.size [%] != cols [%]",
				coeffs.size,
				cols
			).throw
		};

		// NOTE: .asList added to force Collection:flop.
		// Array:flop uses a primitive that has a GC bug:
		// https://github.com/supercollider/supercollider/issues/3454
		// --- fix has been merged, but not in public distro yet.
		// --- TODO: a test of this fix can be tried in the local ATK repo file index_test.scd
		^cols.collect({ |i|
			this.colAt(i) * coeffs[i]
		}).asList.flop.collect(_.sum);
	}


	printOn { | stream |
		if (stream.atLimit) { ^this };
		stream << this.class.name << "[ " ;
		this.printItemsOn(stream);
		stream << " ]" ;
	}

	printItemsOn { | stream |
		this.matrix.do { | item, i |
			if (stream.atLimit) { ^this };
			if (i != 0) { stream.comma };
			stream << "\n";
			stream.tab;
			item.printOn(stream);
		};
		stream << "\n";
	}

}