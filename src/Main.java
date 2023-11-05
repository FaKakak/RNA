import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

public class Main {
    public static final int STEPNUM = 150000;
    public static final int INOCUSTEP = 10000;
    public static final int STAREC = 0;
    public static final int RECINT = 5000;

    public static void main(String[] args) {
        long startTime = System.currentTimeMillis();
        // 相当于int main(void)
        int h=0; // 从0开始
        Environment e = new Environment(); // 等价于inits();初始化环境
        // 这个环境里有 int[][] raw_arr; 每格都撒了100个左右的raw_material
        //    RNA[][][] cell_head; 每格都初始化了一个RNA，当作链表头
        Unit_case_generator unitCaseGenerator = new Unit_case_generator(e);
        Record_generator recordGenerator = new Record_generator(e);
        Inoculate_generator inoculateGenerator = new Inoculate_generator(e);
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
        try (PrintWriter writer = new PrintWriter(new FileWriter("pictureJava.txt", true))) {
            writer.printf("程序运行时间："+duration*0.001+" s\n");
        } catch (IOException E) {
            System.out.println("cannot open file");
            System.exit(-1);
        }
    }
}