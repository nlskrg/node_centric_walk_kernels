package datastructure;

public class MatrixTool {
	
	public static <T extends Number> double getElementSum(Matrix<T> m) {
		
		double d = 0;
		
		for (int i=0; i<m.getRowDimension(); i++) {
			for (int j=0; j<m.getColumnDimension(); j++) {
				d += m.get(i, j).doubleValue();
			}
		}
		
		return d;
	}
	
	public static void add(Matrix<Integer> m, Matrix<Integer> a) {
		for (int i=0; i<m.getRowDimension(); i++) {
			for (int j=0; j<m.getColumnDimension(); j++) {
				m.set(i, j, m.get(i, j)+a.get(i, j));
			}
		}
	}
	
	public static void mult(Matrix<Integer> m, int a) {
		for (int i=0; i<m.getRowDimension(); i++) {
			for (int j=0; j<m.getColumnDimension(); j++) {
				m.set(i, j, m.get(i, j)*a);
			}
		}
	}
	
	public static void mult(Matrix<Double> m, double a) {
		for (int i=0; i<m.getRowDimension(); i++) {
			for (int j=0; j<m.getColumnDimension(); j++) {
				m.set(i, j, m.get(i, j)*a);
			}
		}
	}
	
	
	public static Matrix<Long> truncate(Matrix<Long> m) {

		int maxRow = m.getRowDimension();
		for (int i=m.getRowDimension()-1; i>=0; i--) {
			boolean allZero = true;
			for (int j=0; j<m.getColumnDimension(); j++) {
				if (m.get(i, j) != 0) {
					allZero = false;
					break;
				}
			}
			if (allZero) {
				maxRow--;
			} else {
				break;
			}
		}
		
		int maxCol = m.getColumnDimension();
		for (int j=m.getColumnDimension()-1; j>=0; j--) {
			boolean allZero = true;
			for (int i=0; i<m.getRowDimension(); i++) {
				if (m.get(i, j) != 0) {
					allZero = false;
					break;
				}
			}
			if (allZero) {
				maxCol--;
			} else {
				break;
			}
		}
		
		Matrix<Long> tM = new FullMatrix<Long>(maxRow, maxCol);
		for (int i=0; i<maxRow; i++) {
			for (int j=0; j<maxCol; j++) {
				tM.set(i, j, m.get(i, j));
			}
		}
		
		return tM;
	}

}
