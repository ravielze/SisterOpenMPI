## Penjelasan Program

Garis besar kerja program adalah master node akan menerima inputan terlebih dahulu (kernel matriks dan input matriks). Kemudian program berjalan denan urutan pengerjaan sesuai yang diberikan spesifikasi tugas, dengan no 2 dan 3 digabung:

1. Master node akan membagi input matriks ke masing-masing node sejumlah jumlah_input_matriks/jumlah_node dengan openMPI. Pengiriman dan penerimaan menggunakan MPI_Send. Hal yang dikirim ialah integer yang bernilai banyaknya matrix yang dikirim ke tiap node, serta masing masing isi matrix yang akan diproses tiap node.

2. Setelah masing-masing node mendapatkan "beban kerja"-nya mereka akan mengerjakan pekerjaannya dengan skema paralelisasi threading dengan openMP. Pada step ini kami menggabungkan konvolusi dengan mencari selisih cell terbesar dan terkecil tiap matrix, karena kami coba lebih cepat dari apabila dipisah. Idenya adalah setiap menghitung hasil konvolusi tepat 1 cell pada matrix keluaran, ubah nilai max dan min dari matrix tersebut. Hal inilah yang juga memotivasi kami untuk mengadakan skema paralelisasi di atas proses konvolusi, karena proses ini mengakses cell yang sama sehingga kami nilai tidak akan menghasilkan speedup signifikan bila proses tersebut dipecah. Kami sajikan 3 kode berbeda yang berbeda syntax namun esensi mirip, yakni dilakukan paralelisasi pada for (jmlMatrix _ row _ col) karena jika pada salah satu aspek saja ditakutkan terjadi ketidakseimbangan pembagian beban pada salah satu dari 5 atau 16 thread (semakin besar aspek yang dibagi akan semakin seimbang jika dibagi 5 atau 16). Tak lupa juga untuk menambahkan reduction max dan min pada array maxArray dan minArray.

3. Setelah masing-masing node menyelesaikan kerjaannya, mereka akan melakukan sorting masing-masing dan melemparkan hasil sorting tersebut ke "partner"-nya. Sorting yang dilakukan pada masing-masing node kami nilai legal karena jika ingin dipararelkan lagi maka sudah bukan sorting dengan openMPI lagi kategorinya, karena harus menggunakan openMP. Nantinya partner-partner ini akan menggabungkan array-nya sendiri dan array yang diterima, menggunakan fungsi merge_array, selanjutnya array hasil gabungan akan dilemparkan ke partner selanjutnya. Penentuan node mana yang melempar/menerima ialah dari banyaknya bit 0 di akhir representasi biner id node ini. Hal ini terus dilakukan hingga semua array terkumpul di node 0 (master).

## Analisa Perbandingan Runtime

Terdapat 5 test case (termasuk TC pada docs tubes) yang diberikan pada spesifikasi tugas. TC 0 sampai 4 memiliki ukuran kasar (estimasi) berturut-turut (2, 3, 3), (5, 33, 36), (15, 100, 100), (100, 666, 50), (24, 5000, 50), degan masing masing triplet menyatakan ukuran kernel, banyak input, dan ukuran matrix input. Cek kompleksitas sederhana, apabila 3 komponen triplet di atas dimisalkan (x,y,z), maka proram ini memiliki kompleksitas O(x*x*y*z*z + y log(y)). Pada eksekusi program di perangkat lokal dengan spesifikasi (
ISI MAS TITO
), kami mendapati bahwa program paralel kami lebih cepat baru pada eksekusi TC2 dengan 3 node dan 5 thread. Sebetulnya hal ini masih dalam ekspektasi kami, karena pada TC yang kecil, cost yang dibutuhkan untuk membagi dan komunikasi secara umum antar node, serta inisiasi thread, tidak sepadan dengan keuntungan yang didapat dari pekerjaan yang dilakukan 2 sampai 4 kali lipat lebih cepat hanya pada sebagian kecil dari proses penghitungan. Speedup yang dirasakan semakin signifikan seiring naiknya kompleksitas dan ukuran TC. Penghitungan singkat kompleksitas program apabila dijalankan dengan n node dan t thread adalah O(xxyzz/nt + ylog(y)/n), tentunya mengabaikan konstanta "cold start" threading dan komunikasi antar thread. Sehingga dapat diperkirakan faktor speedup antara n hingga nt yang semakin jelas terlihat seiring bertambah besarnya orde x,y,z, tentunya bergantung juga pada perbandingan y dengan variabel lain.

## Perbedaan Hasil Program Paralel dan Program Serial

Dari segi hasil, kami mendapati program paralel kami menghasilkan nilai yang sama dengan program serial. Meskipun demikian, kami pernah mengalami perbedaan, tentunya apabila tidak diberikan beberapa batasan pada parameter pragma jika di dalamnya terdapat variabel yang dapat diakses lebih dari satu thread sekaligus. Sebelum solusi yang sekarang diberikan reduction, karena blok kode max/min terdiri dari 2 operasi yakni compare dan assign, bisa saja saat dijalankan 2 thread sekaligus akan terjadi perbedaan hasil akhir yang tertimpa karena tidak atomik. Misalkan var shared = 1, thread 1 sedang menghitung max(shared,3), sementara thread 2 max(shared, 2), ketika mulai, kedua thread akan masuk ke assigning karena 2,3 > 1, namun bila thread 2 selesai lebih lama, shared akan bernilai 2 sehingga salah. Reduksi yang dilakukan program kita hanya pada max/min karena penjumlahan konvolusi tidak kami paralelkan (berjalan pada maximal 1 thread). Proses sorting juga tidak mungkin menghasilkan perbedaan hasil karena setiap saatnya tidak ada variabel yang dapat diakses 2 node sekaligus.

## Variasi Parameter Node dan Thread

| TC  | Node | Thread | Time     |
| --- | ---- | ------ | -------- |
| 1   | 2    | 5      | 0.018965 |
| 1   | 2    | 16     | 0.018403 |
| 1   | 3    | 5      | 0.028222 |
| 1   | 3    | 16     | 0.028862 |
| 1   | 4    | 5      | 0.033333 |
| 1   | 4    | 16     | 0.143552 |
| 2   | 2    | 5      | 0.617713 |
| 2   | 2    | 16     | 0.625250 |
| 2   | 3    | 5      | 0.432420 |
| 2   | 3    | 16     | 0.442991 |
| 2   | 4    | 5      | 0.335820 |
| 2   | 4    | 16     | 0.469933 |
| 3   | 2    | 5      | 0.581831 |
| 3   | 2    | 16     | 0.588090 |
| 3   | 3    | 5      | 1.225156 |
| 3   | 3    | 16     | 0.819680 |
| 3   | 4    | 5      | 0.394746 |
| 3   | 4    | 16     | 0.415618 |
| 4   | 2    | 5      | 8.986783 |
| 4   | 2    | 16     | 8.919745 |
| 4   | 3    | 5      | 6.302404 |
| 4   | 3    | 16     | 6.315424 |
| 4   | 4    | 5      | 8.260329 |
| 4   | 4    | 16     | 7.885514 |

## Analisa Waktu Eksekusi

Berikut akan dibandingkan perbedaan waktu eksekusi program paralel dengan berbagai parameter, tidak dibandingkan dengan serial.

Jika dilihat dari perbedaan node, terlihat tren yakni penambahan node mempengaruhi waktu eksekusi secara signifikan pada TC, terutama yang berbobot cukup besar. Hal ini membuktikan bahwa hipotesis kompleksitas algoritma yang memiliki faktor speedup antara n dan nt, keduanya dipengaruhi n sebagai jumlah node.

Jika dilihat dari perbedaan thread, sebenarnya tidak terlalu terlihat pengaruh signifikan. Perbedaan dengan pengaruh node ini dapat dikatakan mendukung bahwa kompleksitas algoritma suku ylog(y), yang tidak memiliki faktor speedup t, lebih berbobot dari xxyzz, yang dicerminkan dari TC4 dengan jumlah matrix (y) yang dibengkakkan (dari TC123). Hal yang perlu diperhatikan juga ialah speedup dari perbanyakan thread di suatu TC dengan jumlah node n, dapat juga bersifat slowdown bila jumlah node lebih dari n. Ini karena karena beban kerja yang dikerjakan masing-masing node berkurang sehingga cost bottleneck yang dirasakan setiap node saat inisialisasi thread lebih berasa jika task yang dibagi lebih kecil.

Jika dilihat dari perubahan TC, terlihat bahwa terdapat tren dimana TC1 sangat cepat, TC 2 dan 3 serupa namun cenderung lebih cepat pengerjaan 3, serta TC4 yang paling lama. Hal ini, terutama TC2 lebih lama dari TC3, sangat mendukung hipotesis kompleksitas algoritma kami lagi. Ini karena TC2 dan 3 memiliki (x,y,z) = (15, 100, 100), (10, 666, 50). Nilai y yang membengkak pada TC2->3 menyeimbangi x dan z yang turun. Karena derajat x dan z lebih tinggi juga (kuadrat), maka rasio yang rasanya 1:2 dapat mengimbangi pembengkakan y.

## Author

-   13519002 - Steven Nataniel Kodyat
-   13519007 - Muhammad Tito Prakasa
-   13519044 - Kinantan Arya Bagaspati
