#pragma once
#include <list>
#include <thread>
#include <mutex>
#include <condition_variable>

template<class item>
class channel {
private:
	std::list<item> queue;
	std::mutex m;
	std::condition_variable cv;
	std::condition_variable cv2;
	int maxBufferSize;
	bool closed;
public:
	channel(int maxBufferSize) : maxBufferSize(maxBufferSize), closed(false) { }

	void close() {
		std::unique_lock<std::mutex> lock(m);
		closed = true;
		cv.notify_all();
	}
	bool is_closed() {
		std::unique_lock<std::mutex> lock(m);
		return closed;
	}
	int size() {
		std::unique_lock<std::mutex> lock(m);
		return queue.size();
	}
	void put(const item &i) {
		std::unique_lock<std::mutex> lock(m);
		if(closed)
			throw std::logic_error("put to closed channel");
		if (queue.size() >= 8) {
			cv2.wait(lock, [&](){ return queue.size() < 8; });
		}
		queue.push_back(i);
		cv.notify_one();
	}
	bool get(item &out, bool wait = true) {
		std::unique_lock<std::mutex> lock(m);
		if(wait)
			cv.wait(lock, [&](){ return closed || !queue.empty(); });
		if(queue.empty())
			return false;
		out = queue.front();
		queue.pop_front();
		cv2.notify_one();
		return true;
	}
};
