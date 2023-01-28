namespace FastGoat.UserGroup.Characters;

public class Label : ACell
{
    private string Content { get; }

    public Label()
    {
        Content = "";
    }

    public Label(string c)
    {
        Content = c;
    }

    public override string Display() => Content;
}