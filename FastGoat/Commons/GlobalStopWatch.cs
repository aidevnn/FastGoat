using System.Diagnostics;

namespace FastGoat.Commons;

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

    public static void Time(string label, Action action)
    {
        Console.WriteLine($"# {label} Start");
        sw.Restart();
        action();
        sw.Stop();
        Console.WriteLine($"# {label} Time:{sw.ElapsedMilliseconds} ms");
    }
}