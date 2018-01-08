import java.io.*;
import java.util.regex.*;
import java.util.Date;
import java.text.SimpleDateFormat;
//import java.util.Arrays;

public class TeHmm {
    
    static String[] STATES = {
        "None",
        "DNA",
        "LINE",
        "LTR",
        "SINE",
        "Unknown",
        "Simple",
        "Low"
    };
    
    
    
    //static double OBS_P = -19.0843080232;
    static double OBS_P = -0.65;   // TO DO (if too close to zero, log-prob grows to >>0 b/c of state probs; but if too negative, log-prob gets <<0)
    
    static double[] START_P = {
        Math.log(1.0/8.0),
        Math.log(1.0/8.0), 
        Math.log(1.0/8.0), 
        Math.log(1.0/8.0), 
        Math.log(1.0/8.0), 
        Math.log(1.0/8.0), 
        Math.log(1.0/8.0), 
        Math.log(1.0/8.0), 
    };

    // --fragX 10 --stateX 0.5
    
    static double[] STATE_P = {
        -0.15130363404031322, -0.6325197690706305, -0.6713109682893635, -0.5928115354996938, -0.6905254738108977, -0.6222843531998695, -0.6896103289330874, -0.6821051738538277
    };
    
    static double[][] TRANS_P = {
        {-0.00013309587379702173, -10.58628209580692, -11.89556502082009, -11.227453031208821, -13.47654046859043, -10.120265739571217, -11.274001240875041, -10.30895793810978},
        {-8.143393692290912, -0.0002906913954902773, -20, -20, -20, -20, -20, -20},
        {-8.411970315746078, -20, -0.00022221632588212919, -20, -20, -20, -20, -20},
        {-9.308458919600719, -20, -20, -9.065825224015594e-05, -20, -20, -20, -20},
        {-7.863574761333576, -20, -20, -20, -0.00038457086486380865, -20, -20, -20},
        {-7.838549990994923, -20, -20, -20, -20, -0.00039431801148833595, -20, -20},
        {-5.9609047647476645, -20, -20, -20, -20, -20, -0.002580906489994263, -20},
        {-6.1380884292548785, -20, -20, -20, -20, -20, -20, -0.002161380948782637},
    };
    
            
    public static void main(String [ ] args) {
        try {
            double[][] pred = readPredictions(
             args[0]);
            int startPosition = Integer.parseInt(args[1]);
            int[] path = viterbi(pred, STATES, START_P, STATE_P, TRANS_P, OBS_P);
            printBed(path, "scaffold_1", startPosition, 0);
            //System.out.print(Arrays.toString(path));
        } catch (Exception e) {
            e.printStackTrace(System.err);
        }
    }
    
    static void printBed(int[] path, String scaf, int offset, int bgState) {
        System.err.println( "printBed\t" + new SimpleDateFormat("HH:mm:ss").format(new Date()));
        int prevState = path[0];
        int left = 0;
        for (int i = 0; i < path.length; i++) {
            if (path[i] != prevState) {
                if (prevState != bgState) {
                    System.out.println(toBedLine(scaf, offset, left, i, prevState));
                }
                prevState = path[i];
                left = i;
            }
        }
        if (prevState != bgState) {
            System.out.println(toBedLine(scaf, offset, left, path.length-1, prevState));
        }
    }
    
    static String toBedLine(String scaf, int offset, int left, int i, int prevState) {
        //int start = offset+left;
        StringBuffer line = new StringBuffer();
        line = line.append(scaf).append("\t");
        line = line.append(offset+left).append("\t");
        line = line.append(offset+i).append("\t");
        line = line.append(STATES[prevState]);
        return line.toString();
    }
            
    static double[][] readPredictions(String file) throws Exception {
        
        // count lines
        BufferedReader reader = new BufferedReader(new FileReader(file));
        int nPred = -5;
        while (reader.readLine() != null) nPred++;
        reader.close();
        if (nPred < 2) {
            throw new Exception("At least 2 predictions are required!");
        }
        
        double[][] pred = new double[nPred][STATES.length];
        BufferedReader br = new BufferedReader(new FileReader(file));
        for (int i=0; i<5; i++) br.readLine(); // skip header (5 lines)
        String dataLine0 = br.readLine();
        int dataStartIndex = locateDataStart(dataLine0);
        pred[0] = parseLine(dataLine0, dataStartIndex);
        for (int i=1; i<nPred; i++) { 
            pred[i] = parseLine(br.readLine(), dataStartIndex);
            if ( 1 == i%1000000 ) {
                System.err.println( "readPredictions\t" + i + "\t" + new SimpleDateFormat("HH:mm:ss").format(new Date()));
            }
        }
        return pred;
    }
    
    // 700001     1:None     1:None       *0.986,0.002,0.002,0,0,0.008,0,0.002
    //    281     1:None   7:Simple   +   0.087,0,0,0,0,0.005,*0.716,0.192
    //static final String pattern = "([\\d,\\.\\*]+)\\s*$";
    static int locateDataStart(String line) throws Exception {
        String pattern = "\\S+\\s*$";
        Pattern r = Pattern.compile(pattern);
        Matcher m;
        m = r.matcher(line);
        if (! m.find()) {
            throw new Exception("Unable to parse: " + line);
        }
        return line.indexOf(m.group(0));
    }
    
    
    
    static double[] parseLine(String line, int dataStartIndex) throws Exception {
        String data = line.substring(dataStartIndex).replaceAll("\\*", "");
        String[] fields = data.split(",");
        double[] logprobs = new double[fields.length];
        for (int i=0; i<fields.length; i++) {
            logprobs[i] = Math.log( Double.parseDouble(fields[i]) );
        }
        return logprobs;
    }
    
    /*
     * States will simply be represented as an array of Strings. The indices of the array will be used for all calculations
     * (such as for other array indices) and the Strings themselves only used for reporting.
     * @param predictions 2D array of predicted log-probabilities at each base of each state ([nBases][nStates])
     * @param states array of state names ([nStates])
     * @param startP array of start log-probabilities ([nStates])
     * @param stateP array of state log-probabilities ([nStates])
     * @param transP 2D array of transition log-probabilities ([nStates][nStates]) (from, to)
     * @param obsP log-probability of the observation, which is independent of state, and assumed here to be constant (1 / genome size)
     * @returns array representing optimal state path (integers refer to indices of state name array)
     */
    static int[] viterbi(
            double[][] pred, String[] states, double[] startP, double[] stateP, double[][] transP, double obsP)
            throws Exception
    {
        
        if (pred.length < 2) { throw new Exception("Must have at least 2 predictions"); }
        System.err.println( "viterbi\t" + new SimpleDateFormat("HH:mm:ss").format(new Date()));
        
        double[][] V = new double[pred.length][states.length];
        int[][] paths = new int[pred.length][states.length];
 
        // Initialize base cases (t == 0)
        for (int s = 0; s < states.length; s++) {
            V[0][s] =  startP[s] + pred[s][0];
            paths[0][s] = s;
        }
        
        // Run Viterbi for t > 0
        for (int t=1; t < pred.length; t++) {
            double pSMaxMax = Double.NEGATIVE_INFINITY;
            for (int s = 0; s < states.length; s++) {
                double pSMax = Double.NEGATIVE_INFINITY;
                int s0Max = -1;
                for (int s0 = 0; s0 < states.length; s0++) {
                    double pS = V[t-1][s0] + transP[s0][s] + pred[t][s] - stateP[s] + obsP;
                    if (pS > pSMax) {
                        s0Max = s0;
                        pSMax = pS;
                    }
                    if (pS > pSMaxMax) {    // For error checking
                        pSMaxMax = pS;
                    }
                }
                //(prob, state) = max( (V[t-1][s0] + transP[s0][s] + pred[s][t] - stateP[s] + obsP , s0) for s0 in states );  // P(s|o)/P(s)
                V[t][s] = pSMax;
                paths[t][s] = s0Max;
            }
            if (pSMaxMax == Double.NEGATIVE_INFINITY) {
                throw new Exception("Zero probability at t = " + t);
            }
            if ( t%1000000 == 1 ) {
                System.err.println( "viterbi\t" + t + "\t" + new SimpleDateFormat("HH:mm:ss").format(new Date()));
            }
        }
        double maxLastP = Double.NEGATIVE_INFINITY;
        int lastState = -1;
        for (int s = 0; s < states.length; s++) {
            if (V[pred.length-1][s] > maxLastP) {
                maxLastP = V[pred.length-1][s];
                lastState = s;
            }
        }
        if (lastState == -1) {
            throw new Exception("Probability underrun error (probabiliy less than minimum value of double type).");
        }
        return pathFor(lastState, paths); // follow backpointers to fill array
    }
    //int foo = 0;
    
    static int[] pathFor( int lastState, int[][] paths ) throws Exception {
        System.err.println( "pathFor\t" + new SimpleDateFormat("HH:mm:ss").format(new Date()));
        int[] path = new int[paths.length];
        path[paths.length-1] = lastState;
        int t = -2;
        try {
            for ( t = paths.length-2; t >= 0; t-- ) {
                path[t] = paths[t+1][path[t+1]];
            }
        } catch (Exception e) {
            System.err.println(paths.length);
            System.err.println(t);
            throw e;
        }
        return path;
    }
    
}
