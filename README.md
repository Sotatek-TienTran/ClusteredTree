# ClusteredTree
Yêu cầu phiên bản Java 8 hoặc hơn
Source code gồm những file chính như sau: 
 - City.java: file này lưu thứ tự và tọa độ các thành phố
 - Edge.java: file chưa thông tin về 1 cạnh bao gồm đỉnh bắt đầu và đỉnh kết thúc
 - Cluster.java: file lưu trữ thông tin về các cluster, bao gồm các đỉnh trong cluster đó
 - InitializeChromosome.java: File này chứa các hàm khởi tạo cá thể. Trong file này có 1 hàm là hàm 
     public double[][] primRSTForClusteredTree(int num_Vertex, int num_Cluster, double[][] weight_Matrix, int[][] vertex_In_Cluster, Random rnd)
   Hàm này sẽ trả về giá trị là một ma trận cạnh được khởi tạo là một cây khung với các đỉnh có nối với nhau thì là 1, các đỉnh không nối với nhau là 0. 
 - Chromosome.java: file lưu thông tin về mỗi cá thể, gồm ma trận cạnh (lấy từ hàm khởi tạo ở file InitializeChromosome ở trên), factorialCost, 
   factorialRank, ScalarFitness...
 - ChromosomeCmp.java: file chứa những hàm so sánh 2 chromosome theo cost, theo scalarFitness...để phục vụ sắp xếp quần thể
 - Evaluate.java: File này chứa những hàm đánh giá cho 1 cá thể và hàm decode từ bài toán lớn về bài toán nhỏ
 - Crossover.java, Mutation.java: các file chứa các hàm làm nhiệm vụ lai ghép, đột biến
 - ReadwriteFile.java: đọc file
 - MainClass.java: file chứa hàm đánh giá, sắp xếp quần thể, thực hiện thân chương trình chính 
 - DrawLines.java: file chứa các hàm để vẽ đồ thị ra màn hình
