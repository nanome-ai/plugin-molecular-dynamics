$client = new-object System.Net.WebClient
$client.DownloadFile("http://developer.download.nvidia.com/compute/cuda/10.2/Prod/network_installers/cuda_10.2.89_win10_network.exe","C:\tmp\cuda_10.2.89_win10_network.exe")
cuda_10.2.89_win10_network.exe -s nvcc_10.2 Display.Driver