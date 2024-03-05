import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

public class Record_generator {
    public static final int LONG_CHAIN_LEN = 30;
    private int nr_num=0;
    private int ctr1_num=0;
    private int ctr2_num=0;
    private int ctr3_num=0;
    private int total_mat_num=0;
    private int RNA_num=0;
    private int raw_num=0;
    Environment environment;

    public Record_generator(Environment environment) {
        this.environment = environment;
    }

    public void record(int h,int i){
        for(int y = 0; y< Environment.SIDE; y++){
            for (int x = 0; x < Environment.SIDE; x++) {
                raw_num+=environment.raw_arr[y][x];
                for(RNA rna = environment.cell_head[h][y][x].next;rna.length1!=0;rna=rna.next){
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


            try (PrintWriter writer = new PrintWriter(new FileWriter("JavaModuan.txt", true))) {
                writer.printf("step=%d: nr=%d, ctr1=%d, ctr2=%d, ctr3=%d, RNA=%d  (tn=%d, r=%d)\n", i,
                        nr_num, ctr1_num, ctr2_num, ctr3_num, RNA_num, total_mat_num, raw_num);
            } catch (IOException e) {
                System.out.println("cannot open file");
                System.exit(-1);
            }
        //}

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
