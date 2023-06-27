void synchronize(int threads_count)
{
static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t condvar_in = PTHREAD_COND_INITIALIZER;
static pthread_cond_t condvar_out = PTHREAD_COND_INITIALIZER;
static int threads_in = 0;
static int threads_out = 0;

pthread_mutex_lock(&mutex);
//заблокировали mutex для работы с threads_in и threads_out до сигнала
//дальше идет только один процесс, остальные если есть ждут

threads_in++;
if (threads_in >= threads_count)
{
//текущий поток пришел последним
threads_out = 0;
pthread_cond_broadcast(&condvar_in);
//подали остальным сигнал на разблокировку
}
else
{
while (threads_in < threads_count)
pthread_cond_wait(&condvar_in, &mutex);
//разблокировали mutex (следующий может выполнять)
//и ждем сигнала разблокировки
}
//то же самое для подсчета количества вышедших из фукнции
//(так как возможно использование в цикле)
threads_out++;
if (threads_out >= threads_count)
{
threads_in = 0;
pthread_cond_broadcast(&condvar_out);
}
else
{
while (threads_out < threads_count)
pthread_cond_wait(&condvar_out, &mutex);
}

pthread_mutex_unlock(&mutex);
//все одновременно выходят
}
