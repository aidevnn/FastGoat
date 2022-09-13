using System.Diagnostics;

namespace FastGoat;

public static class GlobalStopWatch
{
    static GlobalStopWatch()
    {
        sw = Stopwatch.StartNew();
    }
    static Stopwatch sw { get; }
    public static void Restart() => sw.Restart();
    public static void Stop() => sw.Stop();
    public static void Show(string label) => Console.WriteLine($"# {label} Time:{sw.ElapsedMilliseconds} ms");
}
