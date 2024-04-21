import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Random;

public class Main {
    public static final int STEPNUM = 3000000;
    public static final int INOCUSTEP = 10000;
    public static final int STAREC = 0;
    public static final int RECINT = 10000;
    public static final int SD = 732;

    public static void main(String[] args) {
        int h=0;
        Random random = new Random(SD);
        Environment e = new Environment(random);
        Unit_case_generator unitCaseGenerator = new Unit_case_generator(e,random);
        Record_generator recordGenerator = new Record_generator(e);
        Inoculate_generator inoculateGenerator = new Inoculate_generator(e,random);
        long startTime = System.currentTimeMillis();
        for (int i = 0; i <= STEPNUM; i++) {
            if(i==INOCUSTEP) inoculateGenerator.inoculate(h);
            if(i>=STAREC&&i%RECINT==0){
                recordGenerator.record(h,i);
            }
            unitCaseGenerator.unit_case(h);
            h=h^1;
        }
        long endTime = System.currentTimeMillis();
        long duration = endTime - startTime;
        try (PrintWriter writer = new PrintWriter(new FileWriter("Java.txt", true))) {
            writer.printf("程序运行时间："+duration*0.001+" s, SD:"+SD+"\n");
        } catch (IOException E) {
            System.out.println("cannot open file");
            System.exit(-1);
        }
    }
}