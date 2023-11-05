import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

public class Record_generator {
    private int nr_num=0;
    private int ctr1_num=0;
    private int ctr2_num=0;
    private int ctr3_num=0;
    private int total_mat_num=0;
    private int RNA_num=0;
    private int raw_num=0;
    Environment environment;
    int[][] raw_arr;
    RNA[][][] cell_head;

    public Record_generator(Environment environment) {
        this.environment = environment;
        this.raw_arr = environment.raw_arr;
        this.cell_head = environment.cell_head;
    }

    public void record(int h,int i){
        for(int y = 0; y< Environment.SIDE; y++){
            for (int x = 0; x < Environment.SIDE; x++) {
                raw_num+=raw_arr[y][x];
                for(RNA rna = cell_head[h][y][x].next;rna.length1!=0;rna=rna.next){
                    RNA_num++;
                    total_mat_num+=rna.length1+rna.length2;
                    if(findSeq(RNA.NRSEQ,rna)) nr_num++;
                    if(findSeq(RNA.CTR1SEQ,rna)) ctr1_num++;
                    if(findSeq(RNA.CTR2SEQ,rna)) ctr2_num++;
                    if(findSeq(RNA.CTR3SEQ,rna)) ctr3_num++;
                }
            }
        }
        total_mat_num+=raw_num;
        System.out.printf("step=%d: nr=%d, ctr1=%d, ctr2=%d, ctr3=%d, RNA=%d  (tn=%d, r=%d)\n", i,
                nr_num, ctr1_num, ctr2_num, ctr3_num, RNA_num, total_mat_num, raw_num);

        if(i==Main.STAREC){
            try (PrintWriter writer = new PrintWriter(new FileWriter("pictureJava.txt"))) {
                writer.printf("step=%d: nr=%d, ctr1=%d, ctr2=%d, ctr3=%d, RNA=%d  (tn=%d, r=%d)\n", i,
                        nr_num, ctr1_num, ctr2_num, ctr3_num, RNA_num, total_mat_num, raw_num);
            } catch (IOException e) {
                System.out.println("cannot open file");
                System.exit(-1);
            }
        }
        else {
            try (PrintWriter writer = new PrintWriter(new FileWriter("pictureJava.txt", true))) {
                writer.printf("step=%d: nr=%d, ctr1=%d, ctr2=%d, ctr3=%d, RNA=%d  (tn=%d, r=%d)\n", i,
                        nr_num, ctr1_num, ctr2_num, ctr3_num, RNA_num, total_mat_num, raw_num);
            } catch (IOException e) {
                System.out.println("cannot open file");
                System.exit(-1);
            }
        }

        nr_num=0;
        ctr1_num=0;
        ctr2_num=0;
        ctr3_num=0;
        total_mat_num=0;
        RNA_num=0;
        raw_num=0;
    }

    public boolean findSeq(char[] subseq,RNA rna){
        boolean hasSeq = false;
        if(rna.length1>= subseq.length){
            for (int i = 0; rna.length1- subseq.length>=i; i++) {
                int j;
                for (j = 0; j < subseq.length; j++) {
                    if(rna.information[0][i+j]!=subseq[j]) break;
                }
                if(j == subseq.length){
                    hasSeq = true;
                    break;
                }
            }
        }
        return hasSeq;
    }
}
