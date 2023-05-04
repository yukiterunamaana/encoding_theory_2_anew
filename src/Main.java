import java.util.Arrays;
import java.util.Random;

class Task1
{
    static final int Q = 1024;
    static final int Q_primal = 1033;
    static final int k=5;
    static final int n=7;
    static final int d = n-k+1;
    static final int t = (int) Math.floor(((double)d-1)/2);

    private static int[] trim(int[] arr) {
        int i = 0;
        while (i < arr.length && arr[i] == 0) {
            i++;
        }
        return Arrays.copyOfRange(arr, i, arr.length);
    }
    public static int[] intToBinary(int n) {
        if (n == 0)
            return new int[] { 0 };
        int[] bits = new int[11];
        int pos = 0;
        while (n > 0) {
            bits[pos++] = n % 2;
            n /= 2;
        }
        int[] result = new int[pos];
        for (int i = 0; i < pos; i++) {
            result[i] = bits[pos - i - 1];
        }
        return result;
    }
    public static int binaryToInt(int[] n) {
        int res = 0;
        for (int i = 0; i < n.length; i++)
            res+=Math.pow(2,n.length-1-i)*n[i];
        return res;
    }
    private static int[] plus_minus_bin_(int[] a, int[] b)
    {
        int[] result = new int[Math.max(a.length, b.length)];
        for (int i = 0; i < result.length; i++) {
            int x = (i < a.length) ? a[i] : 0;
            int y = (i < b.length) ? b[i] : 0;
            result[i] = (x + y) % 2;
        }
        return trim(result);
    }
    private static int[] divide_bin_(int[] a, int[] b)
    {
        if (a.length<b.length) //if deg(a)<deg(b)
            return a;
        else
        {
            int diff = a.length-b.length;
            int[] t = new int[diff+1];
            t[0]=1;

            int[] temp = plus_minus_bin_(a, multiply_bin_(t, b));
            return divide_bin_(temp, b);
        }
    }
    private static int[] multiply_bin_(int[] a, int[] b)
    {
        int[] result = new int[a.length+b.length-1];
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < b.length; j++) {
                result[i+j] += a[i] * b[j];
                result[i + j] %= 2;
            }
        }
        return divide_bin_(result,intToBinary(Q_primal));
    }
    private static int[] power_bin_(int[] a, int p)
    {
        int[] result = {1};
        for (int i = 0; i < p; i++) {
            result=multiply_bin_(result,a);
        }
        return result; //divide_bin_(result,intToBinary(Q_primal));
    }
    public static int plus_minus(int a, int b) {return binaryToInt(plus_minus_bin_(intToBinary(a),intToBinary(b)));}
    public static int divide(int a, int b) {return binaryToInt(divide_bin_(intToBinary(a),intToBinary(b)));}
    public static int multiply(int a, int b) {return binaryToInt(multiply_bin_(intToBinary(a),intToBinary(b)));}
    public static int power(int a, int p) {return binaryToInt(power_bin_(intToBinary(a),p));}

//    public static double[] multiplyPolynomials(double[] a, double[] b) {
//        double[] result = new double[a.length + b.length];
//        Arrays.fill(result, 0);
//        for (int i = 0; i < a.length; i++)
//            for (int j = 0; j < b.length; j++)
//                result[i + j] += a[i] * b[j];
//        int i = result.length - 1;
//        while (i >= 0 && result[i] == 0) i--;
//        double[] trimmed = new double[i+1];
//        System.arraycopy(result, 0, trimmed, 0, i+1);
//        return trimmed;
//    }
//    public static double[] addPolynomials(double[] a, double[] b) {
//        double[] result = new double[Math.max(a.length, b.length)];
//        for (int i = 0; i < result.length; i++) {
//            double x = (i < a.length) ? a[i] : 0;
//            double y = (i < b.length) ? b[i] : 0;
//            result[i] = x + y;
//        }
//        return result;
//    }
//    public static double[] lagrangePolynome(int[] x, int[] y)
//    {
//        double[] res = new double[]{0};
//        for (int i=0;i<x.length;i++)
//            if (y[i]!=0)
//            {
//                double[] member=new double[]{1};
//                for (int j=0; j<x.length; j++)
//                    if (i!=j)
//                    {
//                        double[] x_minus_xj = new double[]{1, -1*x[j]};
//                        member=multiplyPolynomials(member,x_minus_xj);
//                        for (int k=0; k<member.length; k++)
//                        {
//                            if (i!=j)
//                                member[k]/=x[i]-x[j];
//                        }
//                    }
//                for (int k=0; k<member.length; k++)
//                    member[k]*=y[i];
//
//                res=addPolynomials(res,member);
//            }
//        return res;
//    }
    public static int[] addPolynomialsGF(int[] a, int[] b)
    {
        int[] result = new int[Math.max(a.length, b.length)];
        for(int i=0; i<result.length; i++){
            int coeffA = (i<a.length) ? a[i] : 0;
            int coeffB = (i<b.length) ? b[i] : 0;
            result[i] = plus_minus(coeffA, coeffB);
        }
        return result;
    }
    public static int[] multiplyPolynomialsGF(int[] a, int[] b)
    {
        int[] result = new int[a.length + b.length - 1];
        for(int i=0; i<a.length; i++)
            for(int j=0; j<b.length; j++)
                result[i+j] = plus_minus(result[i+j], multiply(a[i], b[j]));
        return result;
    }
    public static int[] lagrangePolynomialsGF(int[] x, int[] y)
    {
        int[] res = new int[]{0};
        for (int i = 0; i < x.length; i++)
            if (y[i] != 0) {
                int[] member = new int[]{1};
                for (int j = 0; j < x.length; j++)
                    if (i != j) {
                        int[] x_minus_xj = new int[]{1, Q - x[j]};
                        member = multiplyPolynomialsGF(member, x_minus_xj);
                        for (int k = 0; k < member.length; k++)
                            if (i != j)
                                member[k] = divide(member[k], plus_minus(x[i], x[j]));
                    }
                for (int k = 0; k < member.length; k++)
                    member[k] = multiply(member[k], y[i]);

                res = addPolynomialsGF(res, member);
            }
        return res;
    }
//    static int calcPolynomeValue(int x, int[] p)
//    {
//        int res=0;
//        for (int i=0; i<p.length; i++)
//            res+=p[p.length-1-i]*Math.pow(x,p.length-1-i);
//        return res % Q;
//    }
    public static int calcPolynomeValueGF(int x, int[] p)
    {
        int res=0;
        for (int i=0; i<p.length; i++)
        //res+=p[p.length-1-i]*Math.pow(x,p.length-1-i);
            res=plus_minus(res,multiply(p[p.length-1-i],power(x,p.length-1-i)));
        return res % Q;
    }
    private static void applyErrors(int[] p)
    {
        Random r = new Random();
        for (int i = 0; i < r.nextInt() % t + 1; i++) //случайное кол-во ошибок, но минимум одна
        {
            int j = Math.abs(r.nextInt() % p.length-1);
            p[j] = Math.abs(r.nextInt() % Q);
        }
    }
    public static int[] gaoEncode(int[] m, int[] a)
    {
        int[] c = new int[a.length];
        for (int i=0;i<a.length;i++)
            c[i]=calcPolynomeValueGF(a[i],m);
        System.out.println("Before error addition:\n"+Arrays.toString(c));
        applyErrors(c);
        return c;
    }
    public int[] extendedEuclidean(int[] g0, int[] g1)
    {
        // Initialize temporary polynomials u and v with the values of g0 and g1 respectively
        int[] u = g0;
        int[] v = g1;
        // Initialize coefficients for the polynomials s and t, where su + tv = gcd(g0, g1)
        int[] s = {1};
        int[] t = {0};
        // Continue the division while the length of add(multiply(u,g0), multiply(v,g1)) is greater than the given constant
        int[] us = multiply_bin_(u,s);
        int[] tv = multiply_bin_(t,v);
        int[] gcd = plus_minus_bin_(us, tv);
        while (gcd.length >= Math.floor((n+k)/2))
        {
//            // Compute the quotient and remainder polynomials
//            int[] quotient = divide_bin_(u, v)[0];
//            int[] remainder = divide_bin_(u, v)[1];
//
//            // Update the temporary polynomials u and v with the remainder and divisor respectively
//            u = v;
//            v = remainder;
//
//            // Update the coefficients for the polynomials s and t
//            int[] temp = plus_minus(s, multiply(quotient, t));
//            s = t;
//            t = temp;
        }

        // Return the coefficients for the polynomial s, which represents the inverse of g0 mod g1
        return s;
    }
    public static int[] gaoDecode(int[] a, int[] b)
    {
        int[] res = new int[]{};

        int[] g0 = new int[] {1};
        for (int i=0; i<a.length; i++)
            g0 = multiply_bin_(g0,new int[]{1, -1*a[i]}); //res*=x-a[i]

        int[] g1 = lagrangePolynomialsGF(a, b);

        return res;
    }
}

public class Main {
    public static void main(String[] args) {
        System.out.printf("d = %d, t = %d\n",Task1.d,Task1.t);
        int[] m = {1,5,3,10,2}; //x4 + 5*x3 + 3*x2 + 10*x + 2
        int[] a = {1,2,3,4,5,6,7};
        System.out.println(Arrays.toString(Task1.gaoEncode(m,a)));
    }
}
