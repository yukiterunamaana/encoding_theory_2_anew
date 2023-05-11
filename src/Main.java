import java.util.*;
import static java.lang.Math.abs;

class Commons
{
    static double[] trim(double[] arr) //V
    {
        int c = 0;
        for (int i = arr.length-1; i>=0; i--)
            if (arr[i]==0.0)
                c++;
            else
                break;
        double[] res = Arrays.copyOfRange(arr,0,arr.length - c);
        return res;
    }
    static int[] trim(int[] arr) //V
    {
        int c = 0;
        for (int i = arr.length-1; i>=0; i--)
            if (arr[i]==0)
                c++;
            else
                break;
        int[] res = Arrays.copyOfRange(arr,0,arr.length - c);
        return res;
    }
    public static double[] addPolynomials(double[] a, double[] b) //V
    {
        double[] result = new double[Math.max(a.length, b.length)];
        int m = Math.min(a.length,b.length);
        for (int i=result.length-1; i>=0; i--)
        {
            result[i]+=(i<a.length)?a[i]:0;
            result[i]+=(i<b.length)?b[i]:0;
        }
        return result;
    }
    public static double[] subPolynomials(double[] a, double[] b) //V
    {
        double[] result = new double[Math.max(a.length, b.length)];
        int m = Math.min(a.length,b.length);
        for (int i = 0; i < result.length; i++)
        {
            result[i] += (i < a.length) ? a[i] : 0;
            result[i] -= (i < b.length) ? b[i] : 0;
        }
        return result;
    }
    public static double[] multiplyPolynomials(double[] a, double[] b) //V
    {
        int n = a.length - 1;
        int m = b.length - 1;
        double[] result = new double[n + m + 1];

        for (int i = 0; i <= n; i++) {
            for (int j = 0; j <= m; j++) {
                result[i + j] += a[i] * b[j];
            }
        }
        return result;
    }
    public static double[] dividePolynomials(double[] a, double[] b, double[] t) //V
    {
        int aDeg = a.length - 1;
        int bDeg = b.length - 1;
        if (aDeg < bDeg) return t;

        double[] quotient = new double[aDeg - bDeg + 1];
        quotient[quotient.length-1]=a[aDeg]/b[bDeg];
        t=addPolynomials(t,quotient);
        double[] sub = multiplyPolynomials(b,quotient);
        double[] newA = subPolynomials(a,sub);
        newA = trim(newA);
        System.out.println(Arrays.toString(a) + " - " + Arrays.toString(sub) + " = " + Arrays.toString(newA));
        return dividePolynomials(newA,b,t);
    }

}

class Task1
{
    static final int Q = 8;
    static final int Q_primal = 11;


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
            result[i] = bits[i];
        }
        return result;
    } //V
    public static int binaryToInt(int[] n) {
        int res = 0;
        for (int i = 0; i < n.length; i++)
            res+=Math.pow(2,i)*n[i];
        return res;
    } //V

//    static int add(int a, int b){return (a+b)%Q;}
//
//    static int sub(int a, int b){return (a-b+Q)%Q;}
//
////    static int pow(int a, int p)
////    {
////        int res=1;
////        for (int i=0; i<p; i++)
////
////    }
//    static int div(int a, int b)
//    {
//        return binaryToInt(divide_bin_(intToBinary(a),intToBinary(b)));
//    }
//
//    static int mult(int a, int b){return div((a%Q)*(b%Q),Q_primal);}

    static int[] plus_minus_bin_(int[] a, int[] b) //V
    {
        int[] result = new int[Math.max(a.length, b.length)];
        for (int i = 0; i < result.length; i++) {
            int x = (i < a.length) ? a[i] : 0;
            int y = (i < b.length) ? b[i] : 0;
            result[i] = (x + y) % 2;
        }
        return Commons.trim(result);
    }

    static int[] divide_bin_(int[] a, int[] b) //V
    {
        int aDeg = a.length - 1;
        int bDeg = b.length - 1;
        if (aDeg < bDeg) return a;
        int[] quotient = new int[aDeg - bDeg + 1];
        quotient[quotient.length-1]=1;
//        int[] sub = multiply_bin_(b,quotient);
        int[] sub = temp_mult_bin(b,quotient);
        int[] temp = plus_minus_bin_(a, sub);
        temp=Commons.trim(temp);
        return divide_bin_(temp, b);
    }
    static int[] temp_mult_bin(int[] a, int[] b) //V
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
    static int[] multiply_bin_(int[] a, int[] b) //V
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
        if (binaryToInt(result)<Q)
            return result;
        else return divide_bin_(result,intToBinary(Q_primal));
    }
    static int[] power_bin_(int[] a, int p) //V
    {
        int[] result = {1};
        for (int i = 0; i < p; i++)
            result=multiply_bin_(result, a);
        return result; //divide_bin_(result,intToBinary(Q_primal));
    }
    public static int plus_minus(int a, int b) //V
    {return binaryToInt(plus_minus_bin_(intToBinary(a),intToBinary(b)));}
    public static int divide(int a, int b) //V
    {return binaryToInt(divide_bin_(intToBinary(a),intToBinary(b)));}
    public static int multiply(int a, int b) //V
    {return binaryToInt(multiply_bin_(intToBinary(a),intToBinary(b)));}
    public static int power(int a, int p) //V
    {return binaryToInt(power_bin_(intToBinary(a),p));}
    public static int[] add_subPolynomialsGF(int[] a, int[] b)
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
        return modulePolynomials_GF(Commons.trim(result),intToBinary(Q_primal));
    }




    public static int[][] dividePolynomials_internal(int[] a, int[] b) {
        int m = a.length - 1;
        int n = b.length - 1;
        if (n == 0 || m < n)
            return new int[][]{new int[]{0},a};
            //throw new IllegalArgumentException("Invalid input polynomial(s)");
        int[] q = new int[m - n + 1];
        int[] r = a;
        for (int i = m; i >= n; i--) {
            int coeff = r[i] / b[n];
            q[i - n] = coeff;
            for (int j = n; j >= 0; j--)
                r[i - n + j] -= coeff * b[j];
        }
        q = Commons.trim(q);
        r = Commons.trim(r);

        return new int[][] {q, r};
        //return new int[][] { modulePolynomials(q,intToBinary(Q_primal)), Commons.trim(r)};
    }
    public static int[] dividePolynomials(int[] polyDividend, int[] polyDivisor)
    {return dividePolynomials_internal(polyDividend, polyDivisor)[0];}
    public static int[] modulePolynomials(int[] polyDividend, int[] polyDivisor)
    {return dividePolynomials_internal(polyDividend, polyDivisor)[1];}

    public static int[][] dividePolynomialsGF_internal(int[] a, int[] b)
    {
        int m = a.length - 1;
        int n = b.length - 1;
        if (n == 0 || m < n)
            throw new IllegalArgumentException("Invalid input polynomial(s)");
        int[] q = new int[m - n + 1];
        int[] r = a;
        for (int i = m; i >= n; i--) {
            int coeff = divide(r[i], b[n]);
            q[i - n] = coeff;
            for (int j = n; j >= 0; j--)
                r[i - n + j] = plus_minus(r[i - n + j],multiply(coeff, b[j]));
        }
        return new int[][] {q, Commons.trim(r)};
    }
    public static int[] dividePolynomials_GF(int[] polyDividend, int[] polyDivisor)
    {return dividePolynomialsGF_internal(polyDividend, polyDivisor)[0];}
    public static int[] modulePolynomials_GF(int[] polyDividend, int[] polyDivisor)
    {return dividePolynomialsGF_internal(polyDividend, polyDivisor)[1];}

//    public static int[] multiplyPolynomialsGF(int[] a, int[] b)
//    {
//        int[] result = new int[a.length + b.length - 1];
//        for(int i=0; i<a.length; i++)
//            for(int j=0; j<b.length; j++)
//                result[i+j] = plus_minus(result[i+j], multiply(a[i], b[j]));
//        return Commons.trim(result);
//    }
//
//    public static int[][] dividePolynomials_internal(int[] polyDividend, int[] polyDivisor) {
//        int dividendDegree = polyDividend.length - 1;
//        int divisorDegree = polyDivisor.length - 1;
//
//        // check for divisor being zero or having higher degree than dividend
//        if (divisorDegree == 0 || dividendDegree < divisorDegree) {
//            throw new IllegalArgumentException("Invalid input polynomial(s)");
//        }
//
//        int[] quotient = new int[dividendDegree - divisorDegree + 1];
//        int[] remainder = Arrays.copyOf(polyDividend, polyDividend.length);
//
//        for (int i = dividendDegree; i >= divisorDegree; i--) {
//            int q = remainder[i] / polyDivisor[divisorDegree];
//            quotient[i - divisorDegree] = q;
//
//            for (int j = divisorDegree; j >= 0; j--) {
//                remainder[i - divisorDegree + j] -= q * polyDivisor[j];
//            }
//        }
//
//        // remove leading zeros if any in remainder
//        int remainderDegree = remainder.length - 1;
//        while (remainderDegree > 0 && remainder[remainderDegree] == 0) {
//            remainderDegree--;
//        }
//        int[] finalRemainder = Arrays.copyOfRange(remainder, 0, remainderDegree + 1);
//
//        return new int[][] { quotient, finalRemainder };
//    }
//
//    public static int[] dividePolynomials(int[] polyDividend, int[] polyDivisor)
//    {return dividePolynomials_internal(polyDividend, polyDivisor)[0];}
//    public static int[] modulePolynomials(int[] polyDividend, int[] polyDivisor)
//    {return dividePolynomials_internal(polyDividend, polyDivisor)[1];}
//
//    public static int[] dividePolynomialsGF(int[] a, int[] b, int[] t)
//    {
//        int aDeg = a.length - 1;
//        int bDeg = b.length - 1;
//        if (aDeg < bDeg) return t;
//
//        int[] quotient = new int[aDeg - bDeg + 1];
//        quotient[quotient.length-1]=divide(a[aDeg],b[bDeg]);
//        t = add_subPolynomialsGF(t,quotient);
//        System.out.println(Arrays.toString(t));
//
//        int[] sub = temp_mult_bin(b,quotient);
//        int[] temp = add_subPolynomialsGF(a, sub);
//        temp=Commons.trim(temp);
//        return dividePolynomialsGF(temp, b, t);
//    }

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
                            member[k] = divide(member[k], plus_minus(x[i], x[j]));
                    }
                for (int k = 0; k < member.length; k++)
                    member[k] = multiply(member[k], y[i]);

                res = add_subPolynomialsGF(res, member);
            }
        return res;
    }
    public static int evalPolynomialGF(int x, int[] p)
    {
        int res=0;
        for (int i=0; i<p.length; i++)
            res=plus_minus(res,multiply(p[p.length-1-i],power(x,p.length-1-i)));
        return res % Q;
    }
    private static void applyErrors(int[] p)
    {
        Random r = new Random();
        int j = abs(r.nextInt() % p.length);
        int pj=p[j];
        while(pj==p[j])
            p[j] = abs(r.nextInt() % Q);
    } //V
    public static int[] gaoEncode(int[] polynome, int[] message)
    {
        int[] c = new int[message.length];
        for (int i=0;i<message.length;i++)
            c[i]=evalPolynomialGF(message[i],polynome);
        System.out.println("Before error addition:\n"+Arrays.toString(c));
        applyErrors(c);
        return c;
    }
//    public static int gcd(int a, int b)
//    {
//        while (b!=0)
//        {
//            int c = a % b;
//            a = b;
//            b = c;
//        }
//        return abs(a);
//    }
//    public static int[] gcd_bin(int[] a, int[] b)
//    {
//        while (b.length>0)
//        {
//            int[] c = modulePolynomials(a,b);
//            a = b;
//            b = c;
//        }
//        return a;
//    }

//    public static int[] gcdPoly(int[] a, int[] b) {
//    int aDeg = a.length - 1;
//    int bDeg = b.length - 1;
//
//    // Check if either polynomial is the zero polynomial
//    if (aDeg < 0) return b;
//    if (bDeg < 0) return a;
//
//    // Set the highest degree polynomial to gcd
//    int[] gcd = (aDeg > bDeg) ? a : b;
//    int[] r = (aDeg > bDeg) ? b : a;
//
//    // Divide until the remainder is zero
//    while (r[0] != 0) {
//        int rDeg = r.length - 1;
//        int gcdDeg = gcd.length - 1;
//        int[] q = new int[gcdDeg - rDeg + 1];
//
//        // Divide the highest degree polynomial of gcd by the highest degree polynomial of remainder
//        q[gcdDeg - rDeg] = gcd[gcdDeg] / r[rDeg];
//
//        // Store the quotient in a new array
//        for (int i = rDeg - 1; i >= 0; i--) {
//            int tmp = 0;
//            for (int j = gcdDeg - rDeg + i; j >= gcdDeg - rDeg; j--) {
//                tmp = gcd[j] - q[gcdDeg - rDeg] * r[j - rDeg + i];
//                gcd[j] = tmp;
//            }
//            gcd[gcdDeg - rDeg + i] = tmp;
//        }
//
//        // Swap gcd and remainder
//        int[] tmp = gcd;
//        gcd = r;
//        r = tmp;
//    }
//
//    return gcd;
//}
public static int[] subPolynomials(int[] a, int[] b)
{
    int[] result = new int[Math.max(a.length, b.length)];
    for(int i=0; i<result.length; i++){
        int coeffA = (i<a.length) ? a[i] : 0;
        int coeffB = (i<b.length) ? b[i] : 0;
        result[i] = coeffA-coeffB;
    }
    return result;
}
public static int[] multiplyPolynomials(int[] a, int[] b)
{
    int[] result = new int[a.length + b.length - 1];
    for(int i=0; i<a.length; i++)
        for(int j=0; j<b.length; j++)
            result[i+j] += a[i]*b[j];
    return Commons.trim(result);
}
//    public static int[] gcd(int[] polynomial1, int[] polynomial2) {
//    // Ensure that polynomial1 is the higher degree polynomial
//    if (polynomial1.length < polynomial2.length) {
//        int[] temp = polynomial1;
//        polynomial1 = polynomial2;
//        polynomial2 = temp;
//    }
//
//    // Keep dividing the lower degree polynomial by the remainder until the remainder is zero
//    int[] remainder = polynomial2;
//    while (!Arrays.equals(remainder, new int[remainder.length])) {
//        int[] quotient = new int[polynomial1.length - remainder.length + 1];
//        quotient[quotient.length - 1] = polynomial1[polynomial1.length - 1] / remainder[remainder.length - 1];
//        int[] temp = multiplyPolynomials(quotient, remainder);
//        polynomial1 = subPolynomials(polynomial1, temp);
//
//        // Swap the polynomials so that the one with the lower degree is the remainder for the next iteration
//        int[] temp2 = Commons.trim(polynomial1);
//        polynomial1 = Commons.trim(remainder);
//        remainder = Commons.trim(temp2);
//    }
//
//    // Return the final remainder (which is the GCD)
//    return polynomial1;
//    }
//    public static int[] gcd_poly(int[] a, int[] b)
//    {
//        int[] c = new int[]{};
//        while (b.length>0)
//        {
//            c = modulePolynomials(a,b);
//            a = b;
//            b = c;
//        }
//        for (int i=0; i<c.length; i++)
//            c[i]/=c[c.length-1];
//
//        return c;
//    }
//    public int[] extendedEuclideanGF(int[] g0, int[] g1)
//    {
//        while (g1.length>0)
//        {
//            int[] c = modulePolynomials(g0,g1);
//            g0 = g1;
//            g1 = c;
//        }
//        return g0;
//    }
//    //TODO
////    public int[] extendedEuclideanGF(int[] g0, int[] g1)
////    {
////        // Initialize temporary polynomials u and v with the values of g0 and g1 respectively
////        int[] u = g0;
////        int[] v = g1;
////        // Initialize coefficients for the polynomials s and t, where su + tv = gcd(g0, g1)
////        int[] s = {1};
////        int[] t = {0};
////        // Continue the division while the length of add(multiply(u,g0), multiply(v,g1)) is greater than the given constant
////        int[] us = multiply_bin_(u,s);
////        int[] tv = multiply_bin_(t,v);
////        int[] gcd = plus_minus_bin_(us, tv);
////        while (gcd.length >= Math.floor((g0.length+g1.length)/2))
////        {
//////            // Compute the quotient and remainder polynomials
//////            int[] quotient = divide_bin_(u, v)[0];
//////            int[] remainder = divide_bin_(u, v)[1];
//////
//////            // Update the temporary polynomials u and v with the remainder and divisor respectively
//////            u = v;
//////            v = remainder;
//////
//////            // Update the coefficients for the polynomials s and t
//////            int[] temp = plus_minus(s, multiply(quotient, t));
//////            s = t;
//////            t = temp;
////        }
////
////        // Return the coefficients for the polynomial s, which represents the inverse of g0 mod g1
////        return s;
////    }
//    //TODO


    //TODO
    static int[] partial_gcd_GF(int[] q0, int[] q1, int stop)
    {
        int[] r0 = q0;
        int[] r1=q1;

        int[] s0=Task1.intToBinary(1);
        int[] s1=Task1.intToBinary(0);

        int[] t0=Task1.intToBinary(0);
        int[] t1=Task1.intToBinary(1);

        while (r1.length-1>stop)
        {
            int[] q = Task1.dividePolynomials_GF(q0,q1);

            int[] ts=s0;
            int[] tt=t0;
            r0=r1;
            s0=s1;
            t0=t1;
            r1=Task1.modulePolynomials_GF(q0,q1);
            s1 = add_subPolynomialsGF(ts, multiplyPolynomialsGF(q, s0));
            t1 = add_subPolynomialsGF(tt, multiplyPolynomialsGF(q, t0));
        }
        return r1;
    }

    //TODO
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

class Task2 extends Task1
{
    final static int BASE = 64;
    final static int l = 3;
    static int[] m = new int[l];
    static Map<Integer, Integer> s = new HashMap<Integer, Integer>();

    static void encode(String fingerprint)
    {
        Random r = new Random();
        for (int i=0; i<BASE; i++)
            if (fingerprint.indexOf(i)!=-1)
                s.put(i, evalPolynomialGF(i,m));
            else s.put(i, r.nextInt() % BASE);
    }
}


public class Main {
    public static void main(String[] args) {
//        final int k=5;
//        final int n=7;
//        final int d = n-k+1;
//        final int t = (int) Math.floor(((double)d-1)/2);
//
//        int[] polynome = new int[]{1,5,3,2,4};
//        int[] message = new int[]{1,2,3,4,5,6,7};
//        int[] code = Task1.gaoEncode(polynome,message);
//        System.out.println(Arrays.toString(code));

//        int[] a = {1, 2, 3}; // 3x^2 + 2x + 1
//        int[] b = {5, 4}; // 4x + 5
//        int[] sum = GaloisFieldGF8.add(a, b); // sum = 3x^2 + 6x + 6
//        int[] difference = GaloisFieldGF8.subtract(a, b); // difference = 3x^2 + 6x + 2
//        System.out.println(Arrays.toString(sum));
//        System.out.println(Arrays.toString(difference));
//        int[] product = GaloisFieldGF8.multiply_bin_(a, b); // product = 5x^3 + 18x^2 + 17x + 20
//        System.out.println(Arrays.toString(product));
//        int[] quotient = GaloisFieldGF8.divide_bin_(a, b); // quotient = 2x + 1
//        System.out.println(Arrays.toString(quotient));
//        //int[] remainder = GaloisFieldGF8.modulo(a, b); // remainder = 4x + 3
//        //System.out.println(remainder);
//        int[] power = GaloisFieldGF8.power_bin_(a, 3); // power = 1x^6 + 4x^5 + 2x^4 + 3x^3 + 2x^2 + 5x + 1
//        System.out.println(power);

        System.out.println(Task1.plus_minus(6,3));
        System.out.println(Task1.multiply(6,3));
        System.out.println(Task1.divide(3,6));
        System.out.println(Task1.power(3,2));

    }
}
