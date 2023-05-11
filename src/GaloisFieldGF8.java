public class GaloisFieldGF8 {

    final static int Q=8;
    final static int Qprime=11;
    // addition of two polynomials in GF(8)

//    public static int add (int a, int b)
//    {
//        return (a+b)%Q;
//    }
//    public static int subtract (int a, int b)
//    {
//        return (a-b+Q)%Q;
//    }
//
//    public static int multiply (int a, int b)
//    {
//        return ()
//    }

    public static int[] add(int[] a, int[] b) {
        int n = Math.max(a.length, b.length);
        int[] result = new int[n];
        for (int i = 0; i < n; i++) {
            int ai = i < a.length ? a[i] : 0;
            int bi = i < b.length ? b[i] : 0;
            result[i] = (ai + bi) % Q;
        }
        return result;
    }

    // subtraction of two polynomials in GF(8)
    public static int[] subtract(int[] a, int[] b) {
        int n = Math.max(a.length, b.length);
        int[] result = new int[n];
        for (int i = 0; i < n; i++) {
            int ai = i < a.length ? a[i] : 0;
            int bi = i < b.length ? b[i] : 0;
            result[i] = (ai - bi);
            if (result[i]<0)
                result[i]+=Q;
        }
        return result;
    }


    public static int[] divide_bin_(int[] a, int[] b) //V
    {
        int aDeg = a.length - 1;
        int bDeg = b.length - 1;
        if (aDeg < bDeg) return a;
        int[] quotient = new int[aDeg - bDeg + 1];
        quotient[quotient.length-1]=a[aDeg]/b[bDeg]; //FIX
//        int[] sub = multiply_bin_(b,quotient);
        int[] sub = temp_mult_bin(b,quotient);
        int[] temp = subtract(a, sub);
        temp=Commons.trim(temp);
        return divide_bin_(temp, b);
    }
    public static int[] temp_mult_bin(int[] a, int[] b) //V
    {
        int n = a.length - 1;
        int m = b.length - 1;
        int[] result = new int[n + m + 1];
        for (int i = 0; i <= n; i++) {
            for (int j = 0; j <= m; j++) {
                result[i + j] += a[i] * b[j];
                result[i + j] %= 2;
            }
        }
        return result;
    }
    public static int[] multiply_bin_(int[] a, int[] b) //V
    {
        int n = a.length - 1;
        int m = b.length - 1;
        int[] result = new int[n + m + 1];
        for (int i = 0; i <= n; i++) {
            for (int j = 0; j <= m; j++) {
                result[i + j] += a[i] * b[j];
                result[i + j] %= 2;
            }
        }
        if (Task1.binaryToInt(result)<Q)
            return result;
        else return divide_bin_(result,Task1.intToBinary(Qprime));
    }
    public static int[] power_bin_(int[] a, int p) //V
    {
        int[] result = {1};
        for (int i = 0; i < p; i++)
            result=multiply_bin_(result, a);
        return result; //divide_bin_(result,intToBinary(Q_primal));
    }


    // multiplication of two polynomials in GF(8)
    public static int[] multiply(int[] a, int[] b) {
        int[] result = new int[a.length + b.length - 1];
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < b.length; j++) {
                result[i+j] ^= GaloisFieldGF8.multiply(a[i], b[j]);
            }
        }
        return result;
    }

    // division of two polynomials in GF(8)
    private static int[][] divide_internal(int[] dividend, int[] divisor) {
        int n = divisor.length;
        int[] remainder = dividend;
        int[] quotient = new int[dividend.length - n + 1];
        for (int i = quotient.length - 1; i >= 0; i--) {
            int factor = GaloisFieldGF8.divide(remainder[n + i - 1], divisor[n - 1]);
            quotient[i] = factor;
            for (int j = 0; j < n; j++) {
                remainder[i + j] ^= GaloisFieldGF8.multiply(factor, divisor[j]);
            }
        }
        return new int[][]{quotient,remainder};
    }

    public static int[] divide(int[] dividend, int[] divisor) {
        return divide_internal(dividend,divisor)[0];
    }
    public static int[] modulo(int[] dividend, int[] divisor) {
        return divide_internal(dividend,divisor)[1];
    }

    // power of a polynomial in GF(8)
    public static int[] power(int[] a, int b) {
        int[] result = {1};
        while (b > 0) {
            if ((b & 1) == 1) {
                result = GaloisFieldGF8.multiply(result, a);
            }
            a = GaloisFieldGF8.multiply(a, a);
            b >>= 1;
        }
        return result;
    }

    // multiplication of two elements in GF(8)
    public static int multiply(int a, int b) {
        if (a == 0 || b == 0) {
            return 0;
        }
        int loga = GaloisFieldGF8.log(a%Q);
        int logb = GaloisFieldGF8.log(b%Q);
        int logab = (loga + logb) % Q;
        return GaloisFieldGF8.exp(logab);
    }

    // division of two elements in GF(8)
    public static int divide(int a, int b) {
        if (b == 0) {
            throw new ArithmeticException("Division by zero");
        }
        int loga = GaloisFieldGF8.log(a);
        int logb = GaloisFieldGF8.log(b);
        int logab = (loga - logb + 7) % 7;
        return GaloisFieldGF8.exp(logab);
    }

    // logarithm of an element in GF(8)
    public static int log(int a) {
        if (a == 0) {
            return -1;
        }
        int result = 0;
        while (a != 1) {
            a = GaloisFieldGF8.exp(a);
            result++;
        }
        return result;
    }

    // exponentiation of an element in GF(8)
    public static int exp(int i) {
        return GaloisFieldGF8.exps[i];
    }

    // precomputed exponents in GF(8)
    private static final int[] exps = {1, 2, 4, 3, 6, 7, 5};

}
