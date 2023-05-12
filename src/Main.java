import java.util.*;

class Commons
{
    static int[] trim(int[] arr) //V
    {
        int c = 0;
        for (int i = arr.length-1; i>=0; i--)
            if (arr[i]==0)
                c++;
            else
                break;
        int[] res = Arrays.copyOfRange(arr,0,arr.length - c);
        if (res.length>0)
            return res;
        else return new int[] {0};
    }

}

class Task1
{
    static final int Q = 8;
    static final int Q_primal = 11;

    static final int n = 5;
    static final int k = 7;

    public static int[] intToBinary(int n0) {
        int n=n0;
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
    static int[][] divide_bin_(int[] a, int[] b, int[] q) //V
    {
        int aDeg = a.length - 1;
        int bDeg = b.length - 1;
        if (aDeg < bDeg) return new int[][]{q,a};
        int[] quotient = new int[aDeg - bDeg + 1];
        quotient[quotient.length-1]=1;
        q=plus_minus_bin_(q,quotient);
        int[] sub = temp_mult_bin(b,quotient);
        int[] temp = plus_minus_bin_(a, sub);
        temp=Commons.trim(temp);
        return divide_bin_(temp, b, q);
    }
    static int[] divide_pow_bin_(int[] a, int[] b) //V
    {
        int pow_a=0;
        int pow_b=0;
        int btoi_a=binaryToInt(a);
        int btoi_b=binaryToInt(b);
        ;
        while (power(2,pow_a)!=btoi_a)
            pow_a++;
        while (power(2,pow_b)!=btoi_b)
            pow_b++;
        int pow_res = Q-1+((pow_a-pow_b)%(Q-1));
        return power_bin_(new int[]{0,1},pow_res);
    }
    static int[] div_bin_(int[] a, int[] b)
    {return divide_bin_(a,b,new int[]{0})[0];}
    static int[] mod_bin_(int[] a, int[] b)
    {return divide_bin_(a,b,new int[]{0})[1];}

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
            }
        }

        for (int k=0; k<result.length; k++)
            result[k] %= 2;

        if (binaryToInt(result)<Q)
            return result;
        else return mod_bin_(result,intToBinary(Q_primal));
    }

    static int[] power_bin_(int[] a, int p)
    {
        int[] result = {1};
        for (int i = 0; i < p; i++)
            result=multiply_bin_(result, a);
        if (binaryToInt(result)<Q)
            return result;
        else return mod_bin_(result,intToBinary(Q_primal));
    }

    static int[] mult_pow_bin_(int[] a, int[] b) //V
    {
        int pow_a=0;
        int pow_b=0;
        int btoi_a=binaryToInt(a);
        int btoi_b=binaryToInt(b);
        while (power(2,pow_a)!=btoi_a)
            pow_a++;
        while (power(2,pow_b)!=btoi_b)
            pow_b++;
        int pow_res = ((pow_a+pow_b)%(Q-1));
        return power_bin_(new int[]{0,1},pow_res);
    }


    public static int plus_minus(int a, int b)
    {return binaryToInt(plus_minus_bin_(intToBinary(a),intToBinary(b)));}
    public static int divide(int a, int b)
    {return binaryToInt(div_bin_(intToBinary(a),intToBinary(b)));}
//    public static int modulo(int a, int b)
//    {return binaryToInt(mod_bin_(intToBinary(a),intToBinary(b)));}
    public static int multiply(int a, int b)
    {return binaryToInt(multiply_bin_(intToBinary(a),intToBinary(b)));}
    public static int power(int a, int p)
    {return binaryToInt(power_bin_(intToBinary(a),p));}
    public static int divpowbin(int a, int b)
    {return binaryToInt(divide_pow_bin_(intToBinary(a),intToBinary(b)));}

    public static int[] add_subPolynomialsGF(int[] a, int[] b) //V
    {
        int[] result = new int[Math.max(a.length, b.length)];
        for(int i=0; i<result.length; i++){
            int coeffA = (i<a.length) ? a[i] : 0;
            int coeffB = (i<b.length) ? b[i] : 0;
            result[i] = plus_minus(coeffA, coeffB);
        }
        return result;
    }

    public static int[] multiplyPolynomialsGF(int[] a, int[] b) //V
    {
        int[] result = new int[a.length + b.length - 1];
        for(int i=0; i<a.length; i++)
            for(int j=0; j<b.length; j++)
                result[i+j] = plus_minus(result[i+j], multiply(a[i], b[j]));
        result=Commons.trim(result);

        //if (binaryToInt(result)<Q)
            return result;
        //else
        //return modulePolynomials_GF(result,intToBinary(Q_primal));
    }
    public static int[][] dividePolynomialsGF_internal(int[] a, int[] b) //V
    {
        int m = a.length - 1;
        int n = b.length - 1;
        if (n == 0 || m < n)
            throw new IllegalArgumentException("Invalid input polynomial(s)");
        int[] q = new int[m - n + 1];
        int[] r = a;
        for (int i = m; i >= n; i--) {
            int coeff = divpowbin(r[i], b[n]);
            q[i - n] = coeff;
            for (int j = n; j >= 0; j--)
                r[i - n + j] = plus_minus(r[i - n + j],multiply(coeff, b[j]));
        }
        r = Commons.trim(r);
        return new int[][] {q, r};
    }
    public static int[] dividePolynomials_GF(int[] polyDividend, int[] polyDivisor)
    {return dividePolynomialsGF_internal(polyDividend, polyDivisor)[0];}
    public static int[] modulePolynomials_GF(int[] polyDividend, int[] polyDivisor)
    {return dividePolynomialsGF_internal(polyDividend, polyDivisor)[1];}

    private static int[] lagrangeMember(int[] x, int[] y, int i)
    {
        int[] member = new int[]{1};
        int coeff = 1;
        for (int j = 0; j < x.length; j++)
            if (i != j) {
                int c = Task1.plus_minus(x[i], x[j]);
                //System.out.println(x[i] + " XOR " + x[j] + " = " + c);
                coeff = Task1.multiply(coeff,c);
            }
        //System.out.println("coeff = " + coeff);
        for (int j = 0; j < x.length; j++)
            if (i != j) {
                int[] x_minus_xj = new int[]{Task1.Q-x[j],1};
                //System.out.println(Arrays.toString(member) + " *= " + Arrays.toString(x_minus_xj));
                member = Task1.multiplyPolynomialsGF(member, x_minus_xj);
                //System.out.println(Arrays.toString(member));
            }

        //System.out.println(Arrays.toString(member) + " /= " + coeff);
        for (int k = 0; k < member.length; k++)
            member[k] = Task1.divide(member[k], coeff);

        //System.out.println(Arrays.toString(member) + " *= " + y[i]);
        for (int k = 0; k < member.length; k++)
            member[k] = Task1.multiply(member[k], y[i]);

        //System.out.println("Total: " + Arrays.toString(member));

        return member;
    }
    public static int[] lagrangePolynomialsGF(int[] x, int[] y)
    {
        int[] res = new int[]{0};
        for (int i=0; i<x.length; i++)
            res=Task1.add_subPolynomialsGF(res, lagrangeMember(x,y,i));
        return res;
    }

    public static int evalPolynomialGF(int x, int[] p) //V
    {
        int res=0;
        for (int i=0; i<p.length; i++)
        {
            int pw = power(x,i);
            //System.out.println(x + " ** " + i + " = " + pw);

            int a = multiply(p[i],pw);
            //System.out.println(pw + " * " + p[i] + " = " + a);

            //System.out.println(res + " XOR " + a + " = " + plus_minus(res, a));
            res=plus_minus(res,a);
            //System.out.println();
        }
        if (res<Q)
            return res;
        else return binaryToInt(mod_bin_(intToBinary(res),intToBinary(Q_primal)));
    }
    private static void applyErrors(int[] p)
    {
        Random r = new Random();
        int j = Math.abs(r.nextInt() % p.length);
        int pj=p[j];
        while(pj==p[j])
            p[j] = Math.abs(r.nextInt() % Q);
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

    public static int[] partial_gcd_GF_simpl(int[] a, int[] b, double stop)
    {
        if (b.length > a.length) {
            int[] temp = a;
            a = b;
            b = temp;
        }
        while (a.length>=stop && b.length > 0 && b[b.length-1]!=0) {

            int[] r = modulePolynomials_GF(a, b);
            a = b;
            b = r;

            if (b[b.length-1]!=0){
                int c = b[b.length-1];
                for (int i = 0; i < b.length; i++) {
                    b[i] = divide(b[i],c);
                }
            }
        }
        return a;
    }
    //TODO
    static int[][] partial_gcd_GF(int[] g0, int[] g1, double stop)
    {
        int[] u = new int[]{1, 0};
        int[] v = new int[]{0, 1};
        while (g0.length>=stop && g1.length!=0 && g1[g1.length-1]!=0) {
            int[][] quotRem = dividePolynomialsGF_internal(g0, g1);
            int[] q = quotRem[0];
            int []r = quotRem[1];
            int[] tempU = add_subPolynomialsGF(u, multiplyPolynomialsGF(q, v));
            u = v;
            v = tempU;
            g0 = g1;
            g1 = r;
        }
        return new int[][]{g0,v};
    }

    //TODO
    public static int[] gaoDecode(int[] a, int[] b)
    {
        int[] res = new int[]{};
        int[] g0 = new int[] {1};
        for (int i=0; i<a.length; i++)
            g0 = multiplyPolynomialsGF(g0,new int[]{a[i],1}); //res*=x-a[i]

        System.out.println("g0: "+Arrays.toString(g0));
        int[] g1 = Commons.trim(lagrangePolynomialsGF(a, b));
        System.out.println("g1: "+Arrays.toString(g1));
        int[][] gv = partial_gcd_GF(g0,g1,(n+k)/2);
        int[] g = Commons.trim(gv[0]);
        int[] v = Commons.trim(gv[1]);
        System.out.println("g: "+Arrays.toString(g));
        System.out.println("v: "+Arrays.toString(v));

        int[][] fr = dividePolynomialsGF_internal(g,v);
        int[] f = fr[0];
        int[] r = fr[1];
        System.out.println("r: "+Arrays.toString(r));
        System.out.println("f: "+Arrays.toString(f));


        return res;
    }

}

class Task2 extends Task1
{
    final static int BASE = 1024;
    static int[] fingerprint = new int[]{};
    static int[] m = new int[]{};
    static Map<Integer, Integer> s = new HashMap<Integer, Integer>();

    static void vault_encode(int[] fingerprint, int[] polynome)
    {
        m = polynome;
        Random r = new Random();
        for (int i=0; i<BASE; i++)
            if (Arrays.asList(fingerprint).contains(i))
                s.put(i, evalPolynomialGF(i,m));
            else s.put(i, r.nextInt() % BASE);
    }
//
//    static int vault_decode(int[] fingerprint)
//    {
//        int[] key = Task1.lagrangePolynomialsGF()
//    }

}

class Task3 extends Task1
{

    public static String reedMullerEncode(String binaryString) {
        // Define parameters
        final int m = 4;
        final int n = (int)Math.pow(2, m);
        final int r = 2;
        final int k = (int)Math.pow(2, r);

        // Convert binary string to boolean array
        boolean[] input = new boolean[binaryString.length()];
        for (int i = 0; i < binaryString.length(); i++) {
            input[i] = binaryString.charAt(i) == '1';
        }

        // Perform Reed-Muller encoding
        boolean[] output = new boolean[n];
        for (int i = 0; i < n; i++) {
            output[i] = false;
            for (int j = 0; j < k; j++) {
                int idx = i & j;
                boolean term = true;
                for (int l = 0; l < r; l++) {
                    term &= ((idx >> l) & 1) == 1;
                }
                output[i] ^= (term ? input[j] : false);
            }
        }

        // Convert boolean array to binary string
        StringBuilder sb = new StringBuilder();
        for (boolean bit : output) {
            sb.append(bit ? '1' : '0');
        }
        return sb.toString();
    }

}

public class Main {
    public static void main(String[] args) {

        int[] code = new int[] {1,3,4,2,5};

        int[] message = new int[] {0,1,2,3,4};

        //System.out.println(Task1.lagrangePolynomialsGF(code,message));

        int[] encrypted = Task1.gaoEncode(code, message);

        System.out.println(Arrays.toString(encrypted));

        System.out.println(Arrays.toString(Task1.gaoDecode(code,encrypted)));


        System.out.println(Task3.reedMullerEncode("00001111001"));
    }
}
