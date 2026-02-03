import struct

def read_vector_from_file(filename):
    """
    Read C++ std::vector<double> from binary file
    :param filename: filename of .bin file
    :return: list of data
    """
    try:
        with open(filename, 'rb') as inFile:
            # 1. 读取 vector 的大小 (size_t)
            #    'Q' 代表 unsigned long long，它与 C++ 的 size_t 在大多数 64-bit 系统上兼容
            size_bytes = inFile.read(struct.calcsize('Q'))
            if not size_bytes:
                return None
            size = struct.unpack('Q', size_bytes)[0]

            # 2. 读取所有的 double 数据
            #    'd' 代表 C++ 的 double 类型
            #    我们需要读取 'size' 个 'd'
            data_bytes = inFile.read(struct.calcsize('d') * size)
            if len(data_bytes) != struct.calcsize('d') * size:
                print("Warning: File is corrupted or incomplete.")
                return None

            # 3. 将字节流解包成一个 Python tuple
            data_tuple = struct.unpack('d' * size, data_bytes)

            # 4. (可选) 将 tuple 转换为 list，方便后续操作
            return list(data_tuple)

    except FileNotFoundError:
        print(f"Error: File not found: {filename}")
        input("Press ENTER to exit")

    except Exception as e:
        print(f"An error occurred while reading the file: {e}")
        input("Press ENTER to exit")
